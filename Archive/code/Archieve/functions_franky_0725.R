# Define a function for Fourier Transform
fft_compress <- function(x) {
  fft_result <- fft(as.numeric(x))
  return(Re(fft_result))
}

## Define a PAA transformation
paa_transform <- function(series, num_segments = 9) {
  segment_size <- floor(length(series) / num_segments)
  compressed <- numeric(num_segments)
  for (i in 1:num_segments) {
    start_idx <- (i - 1) * segment_size + 1
    end_idx <- min(i * segment_size, length(series))
    compressed[i] <- mean(series[start_idx:end_idx])
  }
  return(compressed)
}

# multiple correlation tests
correlation_tests <- function(data, target_column) {
  test_results <- data.frame(Feature = character(), Correlation = numeric(), P_value = numeric(), stringsAsFactors = FALSE)
  
  for (col in names(data)) {
    if (col != target_column) {
      cor_test <- cor.test(data[[col]], data[[target_column]])
      test_results <- rbind(test_results, data.frame(Feature = col, Correlation = cor_test$estimate, P_value = cor_test$p.value))
    }
  }
  return(test_results)
}

###
PAA <- function(series, segments) {
  n <- length(series)
  segment_size <- n / segments
  paa_series <- numeric(segments)
  
  for (i in 1:segments) {
    start_index <- floor((i - 1) * segment_size) + 1
    end_index <- min(floor(i * segment_size), n)
    paa_series[i] <- mean(series[start_index:end_index])
  }
  
  return(paa_series)
}


perform_pca_regression <- function(data, days_to_use, freq = 5) {
  # Filter and prepare the data
  output_days <- data %>%
    filter(days_since_first <= days_to_use) %>%
    distinct() %>%
    group_by(classified_label) %>%
    slice(-1) %>%
    ungroup() %>%
    select(classified_label, gr, days_since_first)
  
  output_days_wide <- output_days %>%
    pivot_wider(names_from = days_since_first, values_from = gr) %>% 
    left_join(dominant_index, by = "classified_label")
  
  # Convert the data to a matrix for convenience
  data <- as.matrix(output_days_wide[, 2:(days_to_use + 2)])
  target <- as.matrix(output_days_wide[, (days_to_use + 3)])
  
  # Perform PCA
  knots_num <- days_to_use / freq
  pca_result <- prcomp(data, center = TRUE, scale. = TRUE)
  score_matrix <- pca_result$x
  transformed_data <- score_matrix[, 1:knots_num]
  
  loadings_matrix <- pca_result$rotation
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  
  # Prepare the target variable
  target <- target[1:29]
  
  # Create a dataframe for regression
  transformed_df <- as.data.frame(transformed_data)
  transformed_df$target <- target
  
  # Perform linear regression
  model <- lm(target ~ ., data = transformed_df)
  
  return(model)
}

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
  
  # Fit the B-spline without any transformation
  predictor_fd <- Data2fd(scaled_time_points, input, basisobj = bspline_basis)
  
  fitted_matrix <- eval.fd(scaled_time_points, predictor_fd)
  
  # Calculate growth rates
  growth_rates <- matrix(NA, nrow = nrow(fitted_matrix) - 1, ncol = ncol(fitted_matrix))
  for (i in 1:(nrow(fitted_matrix) - 1)) {
    growth_rates[i, ] <- (fitted_matrix[i + 1, ] - fitted_matrix[i, ]) / fitted_matrix[i, ]
  }
  
  return(list(raw = input, fitted = fitted_matrix, gr = growth_rates))
}

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


fit_cubic_spline_and_derivative <- function(df) {
  # Load necessary package
  library(splines)
  
  # Get the number of rows and columns
  n_rows <- nrow(df)
  n_cols <- ncol(df)
  
  # Initialize matrices to store fitted values and derivatives
  fitted_matrix <- matrix(NA, nrow = n_rows, ncol = n_cols)
  derivative_matrix <- matrix(NA, nrow = n_rows, ncol = n_cols)
  
  # Define the sequence of knots (assuming equally spaced knots)
  knots <- quantile(seq(1, n_rows), probs = seq(0, 1, length.out = 10))[-c(1, 10)]
  
  # Loop through each column
  for (i in 1:n_cols) {
    # Extract the column as a time series
    y <- df[, i]
    
    # Fit the cubic B-spline model
    spline_basis <- bs(seq(1, n_rows), knots = knots, degree = 3, intercept = TRUE)
    spline_fit <- lm(y ~ spline_basis)
    
    # Store the fitted values in the matrix
    fitted_matrix[, i] <- fitted(spline_fit)
    
    # Calculate the derivative of the fitted curve
    # Derivative is calculated by multiplying coefficients with the derivative of the basis functions
    basis_derivative <- predict(spline_basis, deriv = 1)
    derivative_matrix[, i] <- basis_derivative %*% coef(spline_fit)[-1]
  }
  
  return(list(fitted = fitted_matrix, gr = derivative_matrix, raw = df))
}


complete_dates <- function(data) {
  data %>%
    dplyr::group_by(
      classified_label
    ) %>%
    tidyr::complete(date = seq.Date(min(date), max(date), by = "day")
    ) %>%
    tidyr::replace_na(
      list(n = 0)
    ) %>%
    dplyr::ungroup()
}
