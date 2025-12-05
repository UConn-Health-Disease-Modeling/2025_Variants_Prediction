# load the packages
library(dplyr)
library(tidyr)
library(mgcv)
library(openxlsx)
library(purrr)
library(fda)
library(forecast)


hosp_url <- "Data/weekly-hospital-admissions-covid-per-million.csv"
daily_est_infec_url <- "Data/daily-new-estimated-infections-of-covid-19.xlsx"

hosp.df           <- read.csv(hosp_url)
daily.estinfec.df <- read.xlsx(daily_est_infec_url)

colnames(hosp.df)       <- c("country", "Code", "time", "hosp.per.million")
colnames(daily.estinfec.df) <- c("country", "Code", "time", "ICL.mean", "IHME.mean", "LSHTM.median", "Y.G.mean", "7.days.avg")

daily.estinfec.df$time <- as.Date(daily.estinfec.df$time, origin = "1899-12-31")

# hosp.df |> filter(country == "South Korea") |> head(10)
# daily.estinfec.df |> filter(country == "South Korea") |> filter(time > "2021-11-13") |> head(10)

daily.estinfec.df <- daily.estinfec.df |> 
  dplyr::select("country", "Code", "time", "7.days.avg")

# using the initial infection patterns to predict the hospitalization stats 
hosp.countries <- hosp.df$country |> unique() # only 35 countries in hospitalization data 

daily.estinfec.df <- daily.estinfec.df |> filter(country %in% hosp.countries)  # 31 countries can be paired
infec.countries <- daily.estinfec.df$country |> unique()
hosp.df <- hosp.df |> filter(country %in% infec.countries)

hosp.df <- hosp.df |> select(-Code)
hosp.df$time <- hosp.df$time |> as.Date()
daily.estinfec.df <- daily.estinfec.df |> select(-Code)

hosp.df <- hosp.df |> 
  group_by(country) |>
  arrange(time) |>
  filter(cumsum(hosp.per.million > 0) > 0) |>
  ungroup()

# calculate the mean hosp rate for each country 
mean.hosp.df <- hosp.df |> 
  group_by(country) |> 
  summarise(mean.hosp.per.million = mean(hosp.per.million)) |> 
  ungroup()

process_days <- function(df, days) {
  df %>%
    group_by(country) %>%
    arrange(time) %>%
    filter(cumsum(`7.days.avg` > 0) > 0) %>%
    slice(1:days) %>%
    mutate(day = row_number()) %>%
    dplyr::select(-time) %>%
    pivot_wider(names_from = day, values_from = `7.days.avg`, names_prefix = "Day_") %>%
    ungroup()
}

estinfec.df.30 <- process_days(daily.estinfec.df, 30)
estinfec.df.45 <- process_days(daily.estinfec.df, 45)
estinfec.df.60 <- process_days(daily.estinfec.df, 60)

fit_and_extract_gam <- function(estinfec_df) {
  fit_gam_with_derivative <- function(row) {
    days <- 1:length(row)
    new_days <- seq(min(days), max(days), length.out = 300)  
    
    gam_model <- gam(row ~ s(days))
    
    fitted_values <- predict(gam_model, newdata = data.frame(days = new_days))
    
    # Use finite differences to approximate the derivative
    deriv_values <- c(NA, diff(fitted_values) / diff(new_days))
    
    list(fitted_values = fitted_values, deriv_values = deriv_values)
  }
  
  estinfec_df.gam <- estinfec_df %>%
    rowwise() %>%
    mutate(fit_result = list(fit_gam_with_derivative(c_across(starts_with("Day_"))))) %>%
    ungroup()
  
  fitted_df <- estinfec_df.gam %>%
    dplyr::select(country, fit_result) %>%
    mutate(fitted_values = map(fit_result, "fitted_values"),
           deriv_values = map(fit_result, "deriv_values")) %>%
    dplyr::select(-fit_result) %>%
    unnest_wider(fitted_values, names_sep = "_") %>%
    unnest_wider(deriv_values, names_sep = "_")
  
  return(fitted_df)
}

# Test with the 30-day data
fitted_df_30 <- fit_and_extract_gam(estinfec.df.30)
fitted_df_45 <- fit_and_extract_gam(estinfec.df.45)
fitted_df_60 <- fit_and_extract_gam(estinfec.df.60)

fitted.infec.30 <- fitted_df_30 |> dplyr::select(country, starts_with("fitted")) |> left_join(mean.hosp.df, by = "country")
gr.infec.30     <- fitted_df_30 |> dplyr::select(country, starts_with("deriv"))  |> left_join(mean.hosp.df, by = "country")
fitted.infec.45 <- fitted_df_45 |> dplyr::select(country, starts_with("fitted")) |> left_join(mean.hosp.df, by = "country")
gr.infec.45     <- fitted_df_45 |> dplyr::select(country, starts_with("deriv"))  |> left_join(mean.hosp.df, by = "country")
fitted.infec.60 <- fitted_df_60 |> dplyr::select(country, starts_with("fitted")) |> left_join(mean.hosp.df, by = "country")
gr.infec.60     <- fitted_df_60 |> dplyr::select(country, starts_with("deriv"))  |> left_join(mean.hosp.df, by = "country")

gr.infec.30[is.na(gr.infec.30)] <- 0
gr.infec.45[is.na(gr.infec.45)] <- 0
gr.infec.60[is.na(gr.infec.60)] <- 0

hosp.rate.list <- list(
  fitted.infec.30 = fitted.infec.30,
  gr.infec.30     = gr.infec.30,
  fitted.infec.45 = fitted.infec.45,
  gr.infec.45     = gr.infec.45,
  fitted.infec.60 = fitted.infec.60,
  gr.infec.60     = gr.infec.60 
)

hosp.rate.url <- "Code/ProcessedData/hosp_rate.rds"
# saveRDS(hosp.rate.list, hosp.rate.url)



hosp.rate.list <- readRDS(hosp.rate.url)
hosp.rate.models <- list()

for (type in names(hosp.rate.list)) {
  
  # type = "fitted.infec.30"
  # type = "gr.infec.30"
  
  data <- hosp.rate.list[[type]]
  
  input <- data[, 2:301] |> t() |> as.matrix()
  output <- data$mean.hosp.per.million
  
  # generate basis
  rangeval <- c(0.1, 30)
  simple_basis <- create.fourier.basis(rangeval,9)
  day.index    <- seq(0.1, by = 0.1, length.out = 300)
  # functional representation
  Smooth       <- smooth.basis(day.index, input, simple_basis)
  Smoothfd <- Smooth$fd
  mdl <- fRegress(output ~ Smoothfd)
  
  hosp.rate.models[[type]] <- mdl
}

# calculate the R sq 
for (type in names(hosp.rate.list)) {
  
  # Extract actual and predicted values
  actual <- hosp.rate.list[[type]]$mean.hosp.per.million
  pred   <- hosp.rate.models[[type]]$yhatfdobj
  
  # Calculate R-squared
  SSE  <- sum((actual - pred)^2)
  SSTO <- sum((actual - mean(actual))^2)
  R_squared <- 1 - (SSE / SSTO)
  
  # Fit a linear model (simple linear regression) between actual and pred
  model <- lm(actual ~ pred)
  
  # Extract the p-value from the model summary
  p_value <- summary(model)$coefficients[2, 4]  # p-value of the slope (pred variable)
  
  # Print R-squared
  cat("The R-squared of modeling", type, "is", R_squared, "\n")
  
  # Print p-value, if less than 0.01, print "< 0.01"
  if (p_value < 0.01) {
    cat("The p-value of modeling", type, "is < 0.01", "\n")
  } else {
    cat("The p-value of modeling", type, "is", p_value, "\n")
  }
  
  cat("\n")
}




################################################################################
# typical time series prediction
full.df <- daily.estinfec.df |> 
  left_join(hosp.df, by = c("country", "time")) |>
  dplyr::select(country, time , `7.days.avg`, hosp.per.million)

full.df[is.na(full.df)] <- 0
colnames(full.df) <- c("country", "time", "days.7.avg", "hosp.per.million")

# Function to calculate R-squared
calculate_r_squared <- function(actual, predicted) {
  SSE <- sum((actual - predicted)^2)
  SSTO <- sum((actual - mean(actual))^2)
  R_squared <- 1 - (SSE / SSTO)
  return(R_squared)
}

results <- data.frame(country = character(), R_squared = numeric(), p = numeric(), q = numeric(), s = numeric(), stringsAsFactors = FALSE)

for (country in unique(full.df$country)) {
  
  # Filter the data for the current country and arrange it by time
  country_data <- full.df %>%
    filter(country == !!country) %>%
    arrange(time)
  
  # Remove rows with NA values in days.7.avg or hosp.per.million
  country_data <- country_data %>% filter(!is.na(days.7.avg) & !is.na(hosp.per.million))
  
  # Extract the response (hosp.per.million) and external regressor (days.7.avg)
  hosp_per_million <- country_data$hosp.per.million
  days_7_avg <- country_data$days.7.avg
  
  # Fit an ARIMA model with days.7.avg as an external regressor
  arima_model <- auto.arima(hosp_per_million, xreg = days_7_avg, max.p = 0, seasonal = TRUE)
  
  # Extract the fitted values (predictions)
  predictions <- fitted(arima_model)
  
  # Calculate R-squared
  R_squared <- calculate_r_squared(hosp_per_million, predictions)
  
  # Extract ARIMA model parameters: p, q, s (seasonal)
  p <- arima_model$arma[1]  # Autoregressive terms (p)
  q <- arima_model$arma[2]  # Moving average terms (q)
  s <- arima_model$arma[5]  # Seasonal period (s)
  
  # Append results to the dataframe
  results <- rbind(results, data.frame(country = country, R_squared = R_squared, p = p, q = q, s = s))
}

results |> arrange(desc(R_squared))




