# set root directory
setwd("/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis")

# source the functions (functions_0706.R)
source('Background/functions_0706.R')

# source my functions
source('Code/functions_franky_0725.R')

library(car)
library(tseries)
library(MASS)
library(forecast)
library(ggplot2)
library(regclass)
library(dplyr)
library(forecast)
library(minpack.lm)  # For nonlinear regression

# load the data
output_url <- "Data/output_SouthKorea2.rds"
output <- readRDS(output_url) %>% 
  group_by(
    classified_label
  ) %>%
  mutate(
    first_date = min(date, na.rm = TRUE)
  ) %>%
  mutate(
    days_since_first = as.integer(date - first_date)
  ) %>%
  dplyr::select(
    -first_date
  ) %>% 
  filter(
    days_since_first != 0
  )


# Create a function to fit and forecast the time series model
fit_forecast <- function(data) {
  ts_data <- ts(data$prevalence, frequency = 1)
  fit <- auto.arima(ts_data)
  forecast(fit, h = 200)
}

data <- output %>% 
  dplyr::select(classified_label, prevalence, days_since_first) %>% 
  dplyr::filter(days_since_first <= 90) %>% 
  dplyr::filter(days_since_first != 0) %>% 
  rename(
    date = days_since_first
  )

data %>% filter(classified_label == "B")

logistic_decay <- function(time, P0, r, K) {
  P0 * exp(-r * time) / (1 + (P0 / K) * (exp(-r * time) - 1))
}

fit_logistic_decay <- function(data) {
  time <- data$date
  prevalence <- data$prevalence
  
  # Fit the logistic decay model with increased decay rate
  start_params <- list(P0 = max(prevalence), r = 0.05, K = 1)  # Increased r value for faster decay
  fit <- nlsLM(prevalence ~ logistic_decay(time, P0, r, K), start = start_params, control = nls.lm.control(maxiter = 1000))
  
  # Forecast for the next 500 days
  future_time <- seq(max(time) + 1, max(time) + 500, by = 1)
  forecasted_prevalence <- predict(fit, newdata = data.frame(time = future_time))
  
  # # Introduce fluctuations by adding random noise
  # set.seed(123)  # For reproducibility
  # noise <- rnorm(length(future_time), mean = 0, sd = 0.05)  # Adjust sd for the desired level of fluctuation
  # forecasted_prevalence <- forecasted_prevalence + noise
  
  # Ensure forecasted prevalence stays within bounds [0, 1]
  forecasted_prevalence <- pmin(pmax(forecasted_prevalence, 0), 1)
  
  data.frame(
    Date = future_time,
    Prevalence = forecasted_prevalence,
    Type = "Forecast"
  )
}

results <- data %>%
  group_by(classified_label) %>%
  do(Forecast = fit_logistic_decay(.))


# Extract and visualize the forecast for a specific classified_label
specific_label <- sample(unique(data$classified_label), 1)
historical_data <- data %>%
  filter(classified_label == specific_label)

forecast_data <- results %>%
  filter(classified_label == specific_label) %>%
  pull(Forecast) %>%
  bind_rows()

# Combine historical data with forecast data
historical_data <- historical_data %>%
  mutate(Date = date, Type = "Historical") %>%
  dplyr::select(Date, Prevalence = prevalence, Type)

combined_data <- bind_rows(historical_data, forecast_data)

# Plot using ggplot2
ggplot(combined_data, aes(x = Date, y = Prevalence, color = Type)) +
  geom_line() +
  labs(title = paste("Forecast and Historical Data for", specific_label),
       x = "Date",
       y = "Prevalence") +
  theme_minimal() +
  scale_color_manual(values = c("Historical" = "black", "Forecast" = "blue"))

