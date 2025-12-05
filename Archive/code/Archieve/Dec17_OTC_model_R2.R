library(dplyr)
library(tidyr)
library(glasso)
library(genlasso)
library(stringr)

# Define file paths
file_paths <- list(
  dominance = "Code/ProcessedData/Nov15_variants_dominance.rds",
  duration  = "Code/ProcessedData/Nov25_variants_duration.rds",
  gr_list   = "Code/ProcessedData/Nov15_variants_gr.rds"
)


dominance <- readRDS(file_paths$dominance)
dominance <- dominance %>% 
  separate(col = country_label, into = c("country", "lineage"), sep = "-")

duration  <- readRDS(file_paths$duration)
gr.list   <- readRDS(file_paths$gr_list)

# Clean up objects starting with "Nov15" or "Nov25"
rm(list = ls(pattern = "^file"))



total_result <- list()

for (element in names(gr.list)) {
  
  
  cat("Processing:", element, "\n")
  
  # --- Process Growth Data ------------------------------------------------------
  
  # Get the growth rate
  # gr.days <- gr.list$days_30
  
  # element = "days_30"
  
  gr.days <- gr.list[[element]]
  
  # Process `gr.days` to extract smoothed sharing growth
  sharing.gr <- lapply(names(gr.days), function(i) {
    data <- gr.days[[i]]
    c(i, data$smoothed_sharing_growth[-1])
  })
  names(sharing.gr) <- names(gr.days)
  
  # Helper function to transform processed list into a data frame
  process_gr_data <- function(gr_list) {
    gr_df <- do.call(rbind, lapply(gr_list, function(x) t(data.frame(x))))
    gr_df <- as.data.frame(gr_df) %>%
      mutate(across(-1, as.numeric))
    colnames(gr_df)[1] <- "country_label"
    return(gr_df)
  }
  
  # Generate final sharing growth data frame
  sharing_gr_df <- process_gr_data(sharing.gr)
  rownames(sharing_gr_df) <- NULL
  
  # --- Compute Over-the-Counter (OTC) Data --------------------------------------
  
  # Define grid span for thresholds
  grid.span <- seq(from = -1, to = 1, length.out = 200)
  
  # Helper function to compute OTC values
  compute_otc <- function(series, grid) {
    sapply(grid, function(x) sum(series > x) / length(series))
  }
  
  # Compute OTC matrix
  sharing_otc <- matrix(nrow = nrow(sharing_gr_df), ncol = length(grid.span) + 1)
  
  for (i in seq_len(nrow(sharing_gr_df))) {
    series <- as.numeric(sharing_gr_df[i, -1])
    sharing_otc[i, 1] <- sharing_gr_df[i, 1]
    sharing_otc[i, 2:ncol(sharing_otc)] <- compute_otc(series, grid.span)
  }
  
  # Convert OTC matrix to data frame and separate country_label
  sharing_otc <- as.data.frame(sharing_otc)
  colnames(sharing_otc) <- c("country_label", paste0("grid_", seq_along(grid.span)))
  
  sharing_otc <- sharing_otc %>%
    separate(country_label, into = c("country", "lineage"), sep = "-")
  
  # --- Combine Dominance Data with OTC Data -------------------------------------
  
  comb.df <- dominance %>%
    left_join(sharing_otc, by = c("country", "lineage"))
  
  comb.df[, 4:ncol(comb.df)] <- lapply(comb.df[, 4:ncol(comb.df)], as.numeric)
  
  # Convert relevant columns to a matrix for model fitting
  comb.share <- as.matrix(comb.df[, 4:ncol(comb.df)])
  
  # --- Perform Fused Lasso Fitting ----------------------------------------------
  
  # Add intercept to the predictor matrix
  predictor_matrix <- cbind(1, comb.share)
  
  # Fit the fused lasso model
  fused_lasso_fit <- fusedlasso1d(
    y = comb.df$mean_share,  # Response variable
    X = predictor_matrix     # Predictor matrix (intercept + OTC data)
  )
  
  # Find the best lambda based on the minimum criterion
  best_lambda_index <- which.min(fused_lasso_fit$lambda)
  best_lambda <- fused_lasso_fit$lambda[best_lambda_index]
  coefficients <- coef(fused_lasso_fit, lambda = best_lambda)$beta
  
  # --- Generate Predictions ----------------------------------------------------
  
  # Calculate predictions
  predictions <- predictor_matrix %*% coefficients
  
  # Combine results into a final data frame
  result <- comb.df[, 1:3] %>%
    cbind(predictions = predictions)
  
  colnames(result)[4] <- "predictions"
  
  # Initialize an output data frame
  output <- result %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(`R^2` = cor(mean_share, predictions)) %>%
    dplyr::ungroup()
  
  # Ensure R^2 is numeric
  output$`R^2` <- as.numeric(output$`R^2`)
  
  total_result[[element]] <- output

}


combined_result <- do.call(cbind, total_result)

list.files("Results/")
write.xlsx(combined_result, "Results/otc.models.xlsx")
