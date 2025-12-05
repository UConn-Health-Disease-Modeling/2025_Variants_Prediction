library(dplyr)
library(tidyr)
library(glasso)
library(genlasso)
library(stringr)



# Load required data files
# This script loads RDS files for variants dominance, duration, and growth rates.

# Define file paths
file_paths <- list(
  dominance = "Code/ProcessedData/Nov15_variants_dominance.rds",
  duration  = "Code/ProcessedData/Nov25_variants_duration.rds",
  gr_list   = "Code/ProcessedData/Nov15_variants_gr.rds"
)

# Load data
dominance <- readRDS(file_paths$dominance)
duration  <- readRDS(file_paths$duration)
gr.list   <- readRDS(file_paths$gr_list)

# Clean up objects starting with "Nov15" or "Nov25"
rm(list = ls(pattern = "^file"))



# get the growth rate
gr.60days <- gr.list$days_32

# Convert to the OTCs
# This script processes growth data from `gr.60days` to create two data frames:
# `sharing_gr_df` for sharing growth and `freq_gr_df` for frequency growth.

# Clone the original list for processing
sharing.gr <- gr.60days
freq.gr <- gr.60days

# Process each element in `gr.60days`
for (i in names(gr.60days)) {
  
  # Extract the data for the current country/region
  data <- gr.60days[[i]]
  
  # Assign processed growth data (excluding the first value)
  sharing.gr[[i]] <- c(i, data$smoothed_sharing_growth[-1])
  freq.gr[[i]] <- c(i, data$smoothed_freq_growth[-1])
}

# Helper function to transform the processed list into a data frame
process_gr_data <- function(gr_list) {
  # Convert the list to a data frame
  gr_df <- do.call(rbind, lapply(gr_list, function(x) t(data.frame(x))))
  gr_df <- as.data.frame(gr_df)
  
  # Convert all columns (except the first one) to numeric
  gr_df <- gr_df %>% mutate(across(-1, as.numeric))
  
  # Rename the first column for consistency
  colnames(gr_df)[1] <- "country_label"
  
  return(gr_df)
}

# Apply the processing to create the final data frames
sharing_gr_df <- process_gr_data(sharing.gr)
# freq_gr_df <- process_gr_data(freq.gr)

# Clean up objects 
rm(list = ls(pattern = "^(gr|.*gr$)"))



# Fit the model on frequency growth rate
# Generate a grid for the span
grid.span <- seq(from = 0, to = 1, length.out = 501)

# Preallocate the matrix for OTC results
sharing_otc <- matrix(nrow = nrow(sharing_gr_df), ncol = length(grid.span))

# Function to compute Occupation Time Curve (OTC) for a series
compute_otc <- function(series, grid) {
  sapply(grid, function(x) sum(series > x) / sum(series > 0))  # Proportion of time values > x
}

# Process each row in the data frame
for (i in seq_len(nrow(sharing_gr_df))) {
  series <- as.numeric(sharing_gr_df[i, -1])
  sharing_otc[i, ] <- compute_otc(series, grid.span)
}

# Set up a 1x2 plotting layout
par(mfrow = c(1, 3))  # Two plots side by side

# Specify the label for the country/variant
selected_label <- "United Kingdom-JN.1"  # You can change this to other labels

# Extract growth rate data for the selected label
growth_rate_data <- sharing_gr_df[dominance$country_label == selected_label, -1] %>% 
  unlist() %>% 
  unname() %>% 
  as.numeric()

# Plot 1: Growth rate over time
plot(
  x = 2:32, 
  y = growth_rate_data, 
  type = "p", 
  col = "blue", 
  lwd = 2, 
  pch = 19,  # Solid circle points
  cex = 0.6, 
  ylab = "Growth Rate", 
  xlab = "Days",
  main = paste("Growth Rate\n", "South Korea-JN.1"),  # Add title
  cex.axis = 0.8,  # Reduce axis text size for better fit
  cex.lab = 0.9, 
  ylim = c(0, 0.2)
)

# Add horizontal red dashed lines
abline(h = 0, col = "red", lty = 2)  # Dashed line at y = 0.02
abline(h = 0.02, col = "red", lty = 2)  # Dashed line at y = 0.02
abline(h = 0.12, col = "red", lty = 2)  # Dashed line at y = 0.12

# Extract OTC data for the selected label
otc_data <- sharing_otc[dominance$country_label == selected_label, ]

# Plot 2: OTC
plot(
  x = grid.span, 
  y = otc_data, 
  type = "l", 
  col = "blue", 
  lwd = 2, 
  ylab = "Proportion of Time", 
  xlab = "Growth Rate",
  main = paste("Occupation Time Curve (OTC)\n", selected_label),
  cex.axis = 0.8,
  cex.lab = 0.9
)



# Reset plotting layout
graphics::par(mfrow = c(1, 1))

# ---- Classify groups based on dominance ----
# Define dominance groups
alpha_duration <- dplyr::filter(duration, grepl("\\(B\\.1\\.1\\.7\\)$", country_label))
alpha_dominance <- dplyr::filter(dominance, grepl("\\(B\\.1\\.1\\.7\\)$", country_label))
alpha_group_otc <- sharing_otc[grepl("\\(B\\.1\\.1\\.7\\)$", dominance$country_label), ]

delta_duration <- dplyr::filter(duration, grepl("\\AY\\.4$", country_label))
delta_dominance <- dplyr::filter(dominance, grepl("\\AY\\.4$", country_label))
delta_group_otc <- sharing_otc[grepl("\\AY\\.4$", dominance$country_label), ]

omicron_BA5_duration <- dplyr::filter(duration, grepl("\\(BA\\.5)$", country_label))
omicron_BA5_dominance <- dplyr::filter(dominance, grepl("\\(BA\\.5)$", country_label))
omicron_BA5_otc <- sharing_otc[grepl("\\(BA\\.5)$", dominance$country_label), ]

omicron_BQ1_duration <- dplyr::filter(duration, grepl("\\BQ\\.1$", country_label))
omicron_BQ1_dominance <- dplyr::filter(dominance, grepl("\\BQ\\.1$", country_label))
omicron_BQ1_otc <- sharing_otc[grepl("\\BQ\\.1$", dominance$country_label), ]

omicron_XBB1_duration <- dplyr::filter(duration, grepl("\\XBB\\.1$", country_label))
omicron_XBB1_dominance <- dplyr::filter(dominance, grepl("\\XBB\\.1$", country_label))
omicron_XBB1_otc <- sharing_otc[grepl("\\XBB\\.1$", dominance$country_label), ]

other_JN1_duration <- dplyr::filter(duration, grepl("\\JN\\.1$", country_label))
other_JN1_dominance <- dplyr::filter(dominance, grepl("\\JN\\.1$", country_label))
other_JN1_otc <- sharing_otc[grepl("\\JN\\.1$", dominance$country_label), ]

# ---- Plot Occupation Time Curves (OTC) ----
# Define x-axis values
x_axis <- grid.span

# Set up colors for different groups
group_colors <- c(
  Alpha = "#1f77b4",  # Blue
  Delta = "#d62728",  # Red
  Omicron_BA5 = "#ff7f0e",  # Orange
  Omicron_BQ1 = "#2ca02c",  # Green
  Omicron_XBB1 = "#9467bd",  # Purple
  JN1 = "#e377c2"  # Pink
)

# # ---- Save Plot as PDF ----
# pdf("Results/OTC_Curves.pdf", width = 10, height = 6)  # Set output file name and dimensions
# 
# # Plot OTC for each group
# graphics::matplot(
#   x_axis, t(alpha_group_otc), type = "l", lty = 1, col = group_colors["Alpha"],
#   xlab = "Growth Rate", ylab = "Probability", main = "Occupation Time Curves (OTC)", lwd = 1
# )
# graphics::matlines(x_axis, t(delta_group_otc), type = "l", lty = 1, col = group_colors["Delta"], lwd = 1)
# graphics::matlines(x_axis, t(omicron_BA5_otc), type = "l", lty = 1, col = group_colors["Omicron_BA5"], lwd = 1)
# graphics::matlines(x_axis, t(omicron_BQ1_otc), type = "l", lty = 1, col = group_colors["Omicron_BQ1"], lwd = 1)
# graphics::matlines(x_axis, t(omicron_XBB1_otc), type = "l", lty = 1, col = group_colors["Omicron_XBB1"], lwd = 1)
# graphics::matlines(x_axis, t(other_JN1_otc), type = "l", lty = 1, col = group_colors["JN1"], lwd = 1)
# 
# # Add legend
# graphics::legend(
#   "topright",
#   legend = names(group_colors),
#   col = unname(group_colors),
#   lty = 1, lwd = 1, cex = 0.8, box.lty = 1
# )
# 
# dev.off()  # Close the PDF device and save the file


# ---- Function to calculate mean for OTC, Duration, and Dominance ----
calculate_means <- function(group_otc, group_duration, group_dominance, label) {
  list(
    otc_mean = colMeans(group_otc, na.rm = TRUE),
    duration_mean = mean(group_duration$duration, na.rm = TRUE),
    dominance_mean = mean(group_dominance$mean_share, na.rm = TRUE),
    label = label
  )
}

# ---- Calculate means for each group ----
alpha_means <- calculate_means(alpha_group_otc, alpha_duration, alpha_dominance, "Alpha")
delta_means <- calculate_means(delta_group_otc, delta_duration, delta_dominance, "Delta")
omicron_BA5_means <- calculate_means(omicron_BA5_otc, omicron_BA5_duration, omicron_BA5_dominance, "Omicron BA.5")
omicron_BQ1_means <- calculate_means(omicron_BQ1_otc, omicron_BQ1_duration, omicron_BQ1_dominance, "Omicron BQ.1")
omicron_XBB1_means <- calculate_means(omicron_XBB1_otc, omicron_XBB1_duration, omicron_XBB1_dominance, "Omicron XBB.1")
other_JN1_means <- calculate_means(other_JN1_otc, other_JN1_duration, other_JN1_dominance, "JN.1")

# ---- Store all group data ----
all_means <- list(
  alpha_means,
  delta_means,
  omicron_BA5_means,
  omicron_BQ1_means,
  omicron_XBB1_means,
  other_JN1_means
)

# ---- Set up colors for different groups ----
group_colors <- c(
  `Alpha(B.1.1.7)` = "#1f77b4",  # Blue
  `Delta(AY.4)` = "#d62728",  # Red
  `Omicron(BA.5)` = "#ff7f0e",  # Orange
  `Omicron(BQ.1)` = "#2ca02c",  # Green
  `Omicron(XBB.1)` = "#9467bd",  # Purple
  `JN.1` = "#e377c2"  # Pink
)

# ---- Generate legend text with dynamic averages ----
all_means[[1]]$label <- "Delta(AY.4)"
all_means[[2]]$label <- "Delta(AY.4)"
all_means[[3]]$label <- "Omicron(BA.5)"
all_means[[4]]$label <- "Omicron(BQ.1)"
all_means[[5]]$label <- "Omicron(XBB.1)"
all_means[[6]]$label <- "JN.1"
legend_text <- sapply(all_means, function(group) {
  sprintf(
    "%s\n Duration: %.0f \n Dominance: %.2f",
    group$label, group$duration_mean, group$dominance_mean
  )
})


# ---- Save Plot as PDF ----
pdf("Manuscripts/plots/Mean_OTC_Curves.pdf", width = 10, height = 6)  # Set output file name and dimensions

# ---- Plot Mean OTC Curves ----
graphics::par(mar = c(5, 4, 4, 12))  # Increase the right margin for the legend
graphics::plot(
  x_axis, alpha_means$otc_mean, type = "l", lty = 1, col = group_colors["Alpha(B.1.1.7)"],
  xlab = "Growth Rate", ylab = "Cumulative Probability", main = "Mean Occupation Time Curves (OTC)", lwd = 2, 
  xlim = c(0, 0.4)
)
graphics::lines(x_axis, delta_means$otc_mean,        type = "l", lty = 1, col = group_colors["Delta(AY.4)"], lwd = 2)
graphics::lines(x_axis, omicron_BA5_means$otc_mean,  type = "l", lty = 1, col = group_colors["Omicron(BA.5)"], lwd = 2)
graphics::lines(x_axis, omicron_BQ1_means$otc_mean,  type = "l", lty = 1, col = group_colors["Omicron(BQ.1)"], lwd = 2)
graphics::lines(x_axis, omicron_XBB1_means$otc_mean, type = "l", lty = 1, col = group_colors["Omicron(XBB.1)"], lwd = 2)
graphics::lines(x_axis, other_JN1_means$otc_mean,    type = "l", lty = 1, col = group_colors["JN.1"], lwd = 2)

# ---- Add Legend ----
graphics::legend(
  "topright",
  inset = c(0, 0),  # Position the legend outside the plot area
  legend = legend_text,
  col = unname(group_colors),
  lty = 1, lwd = 2, cex = 0.5, box.lty = 1, xpd = TRUE,  # xpd = TRUE allows drawing outside the plot area
  y.intersp = 2
)

dev.off()  # Close the PDF device and save the file



# ---- Fit Dominance Model ----

# # Perform fused lasso fitting on dominance data
# fused_lasso_fit <- fusedlasso1d(
#   y = dominance$mean_share,   # Response variable (mean share)
#   X = cbind(1, sharing_otc)   # Predictor matrix (intercept + OTC data)
# )

# # Save the fitted dominance model to an RDS file
# saveRDS(fused_lasso_fit, file = "Results/fused_lasso_fit_mean_share.rds")

fused_lasso_fit <- readRDS("Results/fused_lasso_fit_mean_share.rds")
best_lambda_index <- which.min(fused_lasso_fit$lambda)  # Index of the best lambda
best_lambda <- fused_lasso_fit$lambda[best_lambda_index]  # Value of the best lambda

# Extract coefficients for the best lambda
coefficients <- coef(fused_lasso_fit, lambda = best_lambda)$beta

# ---- Make Predictions ----

# Generate predictions using the fitted model
predictor_matrix <- cbind(1, sharing_otc)  # Add intercept term to the predictors
predictions <- predictor_matrix %*% coefficients  # Calculate predictions

# ---- Output Results ----
# Add predictions to the dominance data
dominance$pred <- predictions

# Define the data source and read the CSV
url <- "UKHSA-UConn-variant-modelling/variant_modelling/data/summary_GISAID_20240918.csv"
data <- read.csv(url)

# Extract main countries (top 30 by total count)
main_countries <- data |>
  group_by(country) |>
  summarise(num_total = sum(numerator)) |>
  arrange(desc(num_total)) |>
  head(30) |>
  dplyr::select(country) |>
  unlist() |>
  as.vector()

# saveRDS(main_countries, "Code/ProcessedData/main_countries.rds")

# Initialize an empty data frame to store results
all_rank_data <- data.frame()

# Loop through each country and extract the top 10 variants
for (country in main_countries) {
  filtered_data <- dominance |>
    dplyr::filter(str_starts(country_label, country)) |>
    dplyr::select(country_label, mean_share, pred) |>
    arrange(desc(mean_share)) 
  
  # Add country column for identification
  filtered_data <- filtered_data |>
    mutate(country = country)
  
  # Combine results into the final data frame
  all_rank_data <- bind_rows(all_rank_data, filtered_data)
}

# Write all top 10 rankings into a single CSV file
all_rank_data[,3][all_rank_data[,3] < 0] <- 0
all_rank_data[,3][all_rank_data[,3] > 1] <- 1
write.csv(all_rank_data, file = "Results/variants_by_country.csv", row.names = FALSE)



# ---- Fit Duration Model ----

# Perform fused lasso fitting for duration
# fused_lasso_fit_duration <- fusedlasso1d(
#   y = duration$duration,    # Response variable (duration)
#   X = cbind(1, sharing_otc) # Predictor matrix (intercept + OTC data)
# )

# # Save the duration dominance model to an RDS file
# saveRDS(fused_lasso_fit, file = "Results/fused_lasso_fit_duration.rds")

# Identify the best lambda value for duration model
fused_lasso_fit <- readRDS("Results/fused_lasso_fit_duration.rds")
best_lambda_index_duration <- which.min(fused_lasso_fit_duration$lambda)
best_lambda_duration <- fused_lasso_fit_duration$lambda[best_lambda_index_duration]

# Extract coefficients for the best lambda
duration_coefficients <- coef(fused_lasso_fit_duration, lambda = best_lambda_duration)$beta

# ---- Load Top 10 Variants Data ----
# Read the previously saved top10_variants_by_country.csv
top10_variants <- read.csv("Results/top10_variants_by_country.csv")

# ---- Prepare Prediction Matrix ----
# Extract the predictor matrix for the country labels in the top 10 data
prediction_matrix <- sharing_otc[match(top10_variants$country_label, dominance$country_label), ]

# Add intercept term
prediction_matrix <- cbind(1, prediction_matrix)

# ---- Make Predictions ----
# Calculate predictions using the fitted duration model
top10_variants$duration_pred <- as.numeric(prediction_matrix %*% duration_coefficients)

# Include the original duration, start, and start_date columns in the output
top10_variants$original_duration <- duration$duration[match(top10_variants$country_label, dominance$country_label)]
top10_variants$start <- duration$start[match(top10_variants$country_label, dominance$country_label)]
top10_variants$start_date <- duration$start_date[match(top10_variants$country_label, dominance$country_label)]

# ---- Save Results ----
# Save predictions and original data into a new CSV file
write.csv(top10_variants, file = "Results/top10_variants_with_duration_predictions.csv", row.names = FALSE)
