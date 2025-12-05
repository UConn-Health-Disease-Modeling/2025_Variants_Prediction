# load the packages
library(dplyr)
library(tidyr)
library(mgcv)
library(openxlsx)
library(purrr)
library(fda)
library(forecast)
library(vars)
library(xgboost)
library(caret)


################################################################################
# processing the infection-hospitalization data (South Korea, United Kingdom, United States)
# assume this data is accurate 

hosp_url <- "Data/weekly-hospital-admissions-covid-per-million.csv"
hosp.df           <- read.csv(hosp_url)

colnames(hosp.df)       <- c("country", "Code", "date", "hosp.per.million")
hosp.df$date <- as.Date(hosp.df$date)

full.df <- hosp.df |> 
  dplyr::select(country, date, hosp.per.million) |> 
  filter(country %in% c("South Korea", "United Kingdom", "United States"))

# hosp.df %>% filter(country == "South Korea")





################################################################################
# processing the variants data 

url <- "UKHSA-UConn-variant-modelling/variant_modelling/data/summary_GISAID_20240918.csv"
data <- read.csv(url)

data2 <- data |> 
  filter(country %in% c("South Korea", "United Kingdom", "USA")) |> 
  mutate(country = ifelse(country == "USA", "United States", country)) |> 
  arrange(country, date)

data3 <- data2 |> 
  rename(n = numerator) |> 
  dplyr::select(country, date, lineage, n)

save_url <- "Code/ProcessedData/global_classification_to_append.rds"
classification_to_append <- readRDS(save_url)

classified_data <- data3 |>
  dplyr::mutate(
    decimal_lineage = paste0(lineage, '.')
  ) |>
  fuzzyjoin::fuzzy_left_join(
    classification_to_append,
    by = 'decimal_lineage',
    match_fun = stringr::str_starts
  ) |>
  dplyr::filter(
    is.na(classified_unasliased) | stringr::str_length(classified_unasliased) == max(stringr::str_length(classified_unasliased)),
    .by = c(lineage, date)
  ) |>
  dplyr::mutate(
    length_class = stringr::str_count(classified_unasliased, '\\.'),
    length_lineage = stringr::str_count(lineage, '\\.'),
    classified_label = dplyr::if_else((length_lineage - length_class) <= 1,
                                      classified_label,
                                      NA),
    classified_label = dplyr::if_else(is.na(classified_label),
                                      "Other",
                                      classified_label)
  ) |>
  dplyr::filter(classified_label != 'Other') |>
  dplyr::select(country, date, lineage, n, classified_label)

# unique(classified_data$classified_label)

classified_data <- classified_data %>% 
  dplyr::select(
    country, date, n, classified_label
  ) %>% 
  dplyr::group_by(
    country, classified_label, date, 
  ) %>% 
  dplyr::summarise(
    n_total = sum (n)
  )

# complete the dates to be continuous
complete_dates <- function(data) {
  data %>%
    dplyr::group_by(country) %>%
    dplyr::mutate(
      min_date = min(date),
      max_date = max(date)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(country, classified_label) %>%
    tidyr::complete(
      date = seq.Date(min(min_date), max(max_date), by = "day")
    ) %>%
    tidyr::replace_na(
      list(n_total = 0)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-min_date, -max_date)  # Drop temporary columns
}

classified_data$date <- as.Date(classified_data$date)

classified_data <- complete_dates(classified_data) %>% 
  dplyr::rename(
    numerator = n_total
  ) %>% 
  dplyr::group_by(
    country, date
  ) %>% 
  dplyr::mutate(
    denominator = sum(numerator, na.rm = TRUE)
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    share = numerator/denominator
  )

(classified_data %>% filter(country == "South Korea"))$classified_label %>% table()


process_variants <- function(data, country_name) {
  variants <- data %>% filter(country == country_name)
  dates <- unique(variants$date)
  variants_list <- lapply(unique(variants$classified_label), function(label) {
    variants %>% filter(classified_label == label) %>% pull(share)
  })
  names(variants_list) <- unique(variants$classified_label)
  variants_df <- cbind(date = dates, as.data.frame(variants_list))
  return(variants_df)
}

# Process each country
variants.South.Korea.df <- process_variants(classified_data, "South Korea")
variants.United.Kingdom.df <- process_variants(classified_data, "United Kingdom")
variants.United.States.df <- process_variants(classified_data, "United States")



full.df.South_Korea <- full.df %>% 
  filter(country == "South Korea") %>% 
  dplyr::select(date, hosp.per.million)
South.Korea.df <- variants.South.Korea.df %>% 
  left_join(full.df.South_Korea, by = "date") %>% 
  filter(!is.na(hosp.per.million))

full.df.United_Kingdom <- full.df %>% 
  filter(country == "United Kingdom") %>% 
  dplyr::select(date, hosp.per.million)
United.Kingdom.df <- variants.United.Kingdom.df %>% 
  left_join(full.df.United_Kingdom, by = "date") %>% 
  filter(!is.na(hosp.per.million))

full.df.United_States <- full.df %>% 
  filter(country == "United States") %>% 
  dplyr::select(date, hosp.per.million)
United.States.df <- variants.United.States.df %>% 
  left_join(full.df.United_States, by = "date") %>% 
  filter(!is.na(hosp.per.million))





rm(list = setdiff(ls(), c("South.Korea.df", "United.Kingdom.df", "United.States.df")))

South.Korea.df <- South.Korea.df[, apply(South.Korea.df, 2, function(col) !all(col == 0))]
United.Kingdom.df <- United.Kingdom.df[, apply(United.Kingdom.df, 2, function(col) !all(col == 0))]
United.States.df <- United.States.df[, apply(United.States.df, 2, function(col) !all(col == 0))]

calculate_zero_percentage <- function(df) {
  zero_percentage <- sapply(df[-1], function(col) {
    mean(col == 0) * 100
  })
  zero_percentage_df <- data.frame(
    Variant = names(zero_percentage),
    Zero_Percentage = zero_percentage
  )
  return(zero_percentage_df)
}

Korea.variants <- calculate_zero_percentage(South.Korea.df) %>% 
  arrange(Zero_Percentage) %>%
  head(11) %>% dplyr::select(Variant) %>% unlist() %>% unname()

UK.variants <- calculate_zero_percentage(United.Kingdom.df) %>% 
  arrange(Zero_Percentage) %>%
  head(11) %>% dplyr::select(Variant) %>% unlist() %>% unname()

US.variants <- calculate_zero_percentage(United.States.df) %>% 
  arrange(Zero_Percentage) %>%
  head(11) %>% dplyr::select(Variant) %>% unlist() %>% unname()



Korea.df <- South.Korea.df %>% dplyr::select(date, Korea.variants)
UK.df    <- United.Kingdom.df %>% dplyr::select(date, UK.variants)
US.df    <- United.States.df %>% dplyr::select(date, US.variants)


colnames(Korea.df)
rm(list = setdiff(ls(), c("Korea.df", "UK.df", "US.df")))

gam_smooth <- function(input_df){
  
  smoothed.df <- input_df
  
  for (col in colnames(input_df)[-1]) {
    
    cat("Processing ", col, "\n")
    
    seq <- input_df[[col]]
    index <- 1:nrow(input_df)
    
    # tune the number of knots in GAM 
    k_values <- seq(5, 100, by = 5)
    
    best_gam <- NULL
    lowest_gcv <- Inf
    
    for (k in k_values) {
      gam_model <- gam(seq ~ s(index, bs = "cs", k = k), method = "REML")
      
      gcv_score <- gam_model$gcv.ubre
      
      if (gcv_score < lowest_gcv) {
        lowest_gcv <- gcv_score
        best_gam <- gam_model
      }
    }
    smoothed.df[[col]] <- predict(best_gam)
  }
  
  return(smoothed.df)
}

smoothed.Korea <- gam_smooth(Korea.df)  
smoothed.UK <- gam_smooth(UK.df)
smoothed.US <- gam_smooth(US.df)

smoothed.Korea$hosp.per.million <- Korea.df$hosp.per.million
smoothed.UK$hosp.per.million <- UK.df$hosp.per.million
smoothed.US$hosp.per.million <- US.df$hosp.per.million



# Define a function to calculate growth rate
calculate_growth_rate <- function(df) {

  growth_rate_df <- df
  
  for (col in colnames(df)[-1]) {
    
    growth_rate_df[[col]] <- c(NA, diff(df[[col]]) / lag(df[[col]], 1)[-1])
  }
  
  return(growth_rate_df)
}
gr.Korea <- calculate_growth_rate(smoothed.Korea)[-1, ]
gr.UK    <- calculate_growth_rate(smoothed.UK)[-1, ]
gr.US    <- calculate_growth_rate(smoothed.US)[-1, ]





################################################################################
# define the model, share growth rate -> hosp growth rate
# using xgboost 

create_lagged_data <- function(gr.df, max_lag){

  seq_length <- dim(gr.df)[1]
  
  output <- list()
  
  for (col in colnames(gr.df)) {
    
    for (i in 0:max_lag) {
      
      name_ <- paste0("lag", i)
      seq_ <- gr.df %>% dplyr::select(all_of(col)) %>% unlist() %>% unname()
      
      output[[col]][[name_]] <- lag(seq_, i)[(1+max_lag):seq_length]
      
    }
    
  }
  
  output.df <- output %>% as.data.frame() %>% 
    dplyr::select(-paste0("date.lag", 1:max_lag))
  output.df$date.lag0 <- as.Date(output.df$date.lag0)

  return(output.df)
}

xgb_lag_tune <- function(df){
  
  model_history <- list()
  R_sq_history  <- list()
  importance_history <- list()
  
  # df <- gr.Korea
  
  for (lag_ in 1:10) {
    
    lag.data <- create_lagged_data(df, lag_)
    
    X <- lag.data %>% 
      dplyr::select(-ends_with("lag0")) %>% 
      dplyr::select(-starts_with("hosp"))%>% as.matrix()
    y <- lag.data$hosp.per.million.lag0
    
    set.seed(123)
    train_index <- createDataPartition(y, p = 0.8, list = FALSE)
    X_train <- X[train_index, ]
    X_test <- X[-train_index, ]
    y_train <- y[train_index]
    y_test <- y[-train_index]
    
    xgb_model <- xgboost(
      data = X_train, 
      label = y_train, 
      nrounds = 100,            # Number of boosting rounds
      objective = "reg:squarederror",  # Regression task
      max_depth = 6,            # Maximum depth of the tree
      eta = 0.3,                # Learning rate
      verbose = 0               # Turn off printing
    )
    
    predictions <- predict(xgb_model, X_test)
    SSE <- sum((y_test - predictions)^2) 
    SST <- sum((y_test - mean(y_test))^2) 
    R_squared <- 1 - (SSE / SST)
    
    model_history[[paste0("lag", lag_)]] <- xgb_model
    R_sq_history[[paste0("lag", lag_)]]  <- R_squared
    importance_history[[paste0("lag", lag_)]] <- xgb.importance(feature_names = colnames(X_train), model = xgb_model)
    
  }

  return(list(model_history = model_history, 
              R_sq_history = R_sq_history, 
              importance_history = importance_history))
  
}




gr.Korea.models <- xgb_lag_tune(gr.Korea)
gr.Korea.models$R_sq_history %>% unlist() %>% as.vector()
xgb.plot.importance(gr.Korea.models$importance_history$lag2)



gr.UK.models <- xgb_lag_tune(gr.UK)
gr.UK.models$R_sq_history %>% unlist() %>% as.vector()
xgb.plot.importance(gr.UK.models$importance_history$lag5)



gr.US.models <- xgb_lag_tune(gr.US)
gr.US.models$R_sq_history %>% unlist() %>% as.vector()
xgb.plot.importance(gr.US.models$importance_history$lag1)




colnames(gr.US)



