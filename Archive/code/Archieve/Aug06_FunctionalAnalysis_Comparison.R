# source the functions (functions_0706.R)
source(glue::glue('/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis/Background/functions_0706.R'))

# load packages 
library(magrittr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(viridis)
library(scales) 
library(fda)

################################################################################
## Classified the data (for South Koera and United Kingdom together)

# get the alias
alias_list <- get_alias()

# set the time range (as large as possible)
months_from_most_recent <- 365 * 6

# load the data and process 
url <- 'Data/summary_gisaid_new.csv'
data <- readr::read_csv(url) |>
  dplyr::filter(
    # country == country_select,
    lineage != 'Unassigned',
    date >= max(date, na.rm = T) - months_from_most_recent
  ) |>
  dplyr::mutate(
    alias = gsub(
      pattern = "[^a-zA-Z]",
      replacement = "",
      x = lineage
    ),
    second_part = gsub(pattern = "[a-zA-Z]",
                       replacement = "",
                       x = lineage)
  ) |>
  dplyr::left_join(alias_list) |>
  dplyr::mutate(
    unaliased_lineage  =
      dplyr::if_else(!is.na(alias_lineage),
                     paste0(alias_lineage, second_part),
                     lineage),
  )

dataset <- data %>% 
  dplyr::select(country, date, lineage, numerator, unaliased_lineage) %>% 
  dplyr::rename(n = numerator)


# # classified the variants
# summarised <- dataset %>%
#   dplyr::group_by(country, unaliased_lineage) %>% 
#   dplyr::summarise(n = sum(n))
# 
# # set p limit
# p_lim <- 0.001
# 
# # set n limit
# n_lim <- 50
# 
# # get the classification to append
# # save the result since the function costs a lot of time
# classification_to_append <- get_classification(
#   lineages = summarised$unaliased_lineage,
#   number_sequences = summarised$n,
#   p_lim = p_lim, n_lim = n_lim, alias_list
# ) |>
#   dplyr::mutate(
#     decimal_lineage = paste0(classified_unasliased, '.')
#   )
#
# saveRDS(classification_to_append, "Data/Aug06_classification_to_append_Combined.rds")

# load the .rds file
append_url <- "Data/Aug06_classification_to_append_Combined.rds"
classification_to_append <- readRDS(append_url)

# append the classification
classified_data <- dataset |>
  dplyr::mutate(
    decimal_lineage = paste0(unaliased_lineage, '.')
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
    length_lineage = stringr::str_count(unaliased_lineage, '\\.'),
    classified_label = dplyr::if_else((length_lineage - length_class) <= 1,
                                      classified_label,
                                      NA),
    classified_label = dplyr::if_else(is.na(classified_label),
                                      "Other",
                                      classified_label)
  ) |>
  dplyr::filter(classified_label != 'Other') |>
  dplyr::select(country, date, lineage, n, classified_label)

# check the number of the classified labels
# length(unique(classified_data$lineage)) # 1608
# length(unique(classified_data$classified_label)) # 88

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
  

head(classified_data)
colnames(classified_data)

# complete the dates to be continuous
complete_dates <- function(data) {
  data %>%
    dplyr::group_by(
      country, 
      classified_label
    ) %>%
    tidyr::complete(date = seq.Date(min(date), max(date), by = "day")
    ) %>%
    tidyr::replace_na(
      list(n_total = 0)
    ) %>%
    dplyr::ungroup()
}

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

# separate data
Korea_data <- classified_data %>%
  dplyr::filter(
    country == "South Korea"
  ) %>%
  dplyr::arrange(
    classified_label, date
  ) %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    dplyr::n() > 30
  ) %>%
  dplyr::mutate(
    sharing = dplyr::coalesce(share, 0),
    gap = as.numeric(date - dplyr::first(date))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-share)

Korea_dominance <- Korea_data %>% 
  dplyr::group_by(
    classified_label
  ) %>% 
  dplyr::summarise(
    dominance = sum(sharing, na.rm = TRUE),
  ) %>% 
  dplyr::arrange(
    dominance
  )

UK_data <- classified_data %>%
  dplyr::filter(
    country == "United Kingdom"
  ) %>%
  dplyr::arrange(
    classified_label, date
  ) %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    dplyr::n() > 30
  ) %>%
  dplyr::mutate(
    sharing = dplyr::coalesce(share, 0),
    gap = as.numeric(date - dplyr::first(date))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-share)

UK_dominance <- UK_data %>% 
  dplyr::group_by(
    classified_label
  ) %>% 
  dplyr::summarise(
    dominance = sum(sharing, na.rm = TRUE),
  ) %>% 
  dplyr::arrange(
    dominance
  )

# length(unique(Korea_data$classified_label)) # 80
# length(unique(UK_data$classified_label)) # 86

# check if the share is good (Bingo!)
# sample_date <- sample(unique(Korea_data$date), 5)
# sample_index <- sample(1:length(sample_date), 1)
# sum((Korea_data %>% dplyr::filter(date == sample_date[sample_index]))$share)

rm(alias_list, classified_data, data, append_url, url)





################################################################################
fda_matrix_prepare <- function(days_input, df, dominance_df) {
  
  result_df <- df %>%
    dplyr::filter(
      gap <= days_input
    ) %>%
    dplyr::group_by(
      classified_label
    ) %>%
    dplyr::filter(
      sum(sharing == 0) <= (days_input-3)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      gap <= days_input
    ) %>%
    dplyr::select(
      classified_label, 
      gap, 
      sharing
    ) %>%
    reshape2::dcast(
      classified_label ~ gap, 
      value.var = "sharing"
    ) %>%
    dplyr::left_join(
      dominance_df, 
      by = "classified_label")
  
  names    <- result_df[, 1] %>% as.matrix()
  input    <- result_df[, 2:(days_input + 2)] %>% as.matrix()
  response <- result_df[, (days_input + 3)] %>% as.matrix()
  
  return(list(names = names, input = input, response = response))
}


#############################
# For Korea data: 
korea_list <- fda_matrix_prepare(days_input = 30, Korea_data, Korea_dominance)
Korea_names <- korea_list$names
Korea_input <- korea_list$input
Korea_response <- korea_list$response[1:length(Korea_names), ]

# For UK data:
uk_list <- fda_matrix_prepare(days_input = 30, UK_data, UK_dominance)
UK_names <- uk_list$names
UK_input <- uk_list$input
UK_response <- uk_list$response[1:length(UK_names), ]

# remove unnecessary terms 
rm(korea_list, uk_list, Korea_data, Korea_dominance, UK_data, UK_dominance)



#############################
BSpline_gr <- function(Region_input) {
  
  num_samples <- nrow(Region_input)
  num_days    <- ncol(Region_input)
  
  time_points <- 0:(num_days-1)
  scaled_time_points <- time_points / max(time_points)  # Scale to [0, 1]
  
  basis_range <- c(0, 1)
  
  input <- t(Region_input)
  
  n_basis <- round(sqrt(nrow(input))) + 3
  degree <- 2  # Cubic B-spline
  
  bspline_basis <- create.bspline.basis(basis_range, n_basis, norder = degree + 1)
  
  # Apply logarithmic transformation to ensure positivity
  log_input <- log(input + 1e-5)  # Add a small constant to avoid log(0)
  
  predictor_fd <- Data2fd(scaled_time_points, log_input, basisobj = bspline_basis)
  
  log_fitted_matrix <- eval.fd(scaled_time_points, predictor_fd)
  
  # Transform back by exponentiating
  fitted_matrix <- exp(log_fitted_matrix)
  
  growth_rates <- matrix(NA, nrow = nrow(fitted_matrix) - 1, ncol = ncol(fitted_matrix))
  for (i in 1:(nrow(fitted_matrix) - 1)) {
    growth_rates[i, ] <- (fitted_matrix[i + 1, ] - fitted_matrix[i, ]) / fitted_matrix[i, ]
  }
  
  return(list(raw = input, fitted = fitted_matrix, gr = growth_rates))
}

# example to show
UK_output <- BSpline_gr(UK_input)
UK_fitted <- UK_output$fitted
UK_gr <- UK_output$gr
UK_raw <- UK_output$raw

# index <- sample(0:50, 1)
# plot(0:30, UK_raw[ , index], type = 'o', col = 'blue', pch = 16, xlab = 'Time', ylab = 'Value',
#      main = "Comparison of Original vs Fitted Values (First Column)")
# lines(0:30, UK_fitted[ , index], col = 'red', lwd = 2)



# Prepare the functional data object
num_days <- nrow(UK_gr)
time_points <- seq(0, 1, length.out = num_days)

n_basis <- 6
basis_range <- c(0, 1)
bspline_basis <- create.bspline.basis(basis_range, n_basis)

UK_fd <- Data2fd(time_points, UK_gr, basisobj = bspline_basis)

# Fit the functional linear model
fRegress_model <- fRegress(UK_response ~ UK_fd)

# Extract the coefficients of the fit (beta functions)
beta_coefs <- fRegress_model$betaestlist[[2]]$fd$coefs

# Create the functional data object for fitted values
fitted_values_fd <- fd(coef = beta_coefs, basisobj = bspline_basis)

# Evaluate fitted values at the original time points
fitted_values <- eval.fd(time_points, fitted_values_fd)

# Display the first few fitted values for verification
print("Fitted values:")
print(head(fitted_values))

# Compare fitted values with actual UK_response values
comparison <- data.frame(
  Time = time_points,
  Actual = UK_response,
  Fitted = fitted_values
)

# Compute the differences
comparison$Difference <- comparison$Actual - comparison$Fitted

# Display the first few rows of comparison
print("Comparison of Actual and Fitted Values:")
print(head(comparison))

# Visualize the comparison
plot(comparison$Time, comparison$Actual, type = 'l', col = 'blue', lwd = 2,
     ylab = 'Response', xlab = 'Time', main = 'Actual vs Fitted Values')
lines(comparison$Time, comparison$Fitted, col = 'red', lwd = 2, lty = 2)
legend("topright", legend = c("Actual", "Fitted"), col = c("blue", "red"),
       lty = 1:2, cex = 0.8)

# Plot the differences
plot(comparison$Time, comparison$Difference, type = 'l', col = 'green', lwd = 2,
     ylab = 'Difference (Actual - Fitted)', xlab = 'Time', main = 'Differences')
abline(h = 0, col = 'gray', lty = 2)









# Sample input data (UK_input)
set.seed(123)  # For reproducibility
UK_input <- matrix(runif(31 * 5, min = 100, max = 200), nrow = 31, ncol = 5)  # Example data

time_points <- 0:30
scaled_time_points <- time_points / max(time_points)  # Scale to [0, 1]

basis_range <- c(0, 1)
input <- t(UK_input)

# Calculate growth rates
growth_rates <- apply(input, 1, function(row) {
  c(NA, diff(row) / row[-length(row)])  # NA for the first point, as it has no previous point
})

# Remove the first row (NA values) from growth_rates for fitting
growth_rates <- growth_rates[-1, ]

# Set the number of basis functions
n_basis <- round(sqrt(nrow(input)))
degree <- 2  # Cubic B-spline

# Create B-spline basis
bspline_basis <- create.bspline.basis(basis_range, n_basis, norder = degree + 1)

# Fit B-spline to growth rates
fd_par <- fdPar(bspline_basis, Lfdobj=2, lambda=0.01)  # Smoothing parameter
growth_rate_fd <- smooth.basis(scaled_time_points[-1], growth_rates, fd_par)$fd

# Plot the B-spline fit for each time series
plot(growth_rate_fd, main="B-Spline Fit to Growth Rates", xlab="Scaled Time", ylab="Growth Rate")



