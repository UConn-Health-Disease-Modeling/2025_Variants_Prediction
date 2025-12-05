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

# # Separate 'country_label' into 'country' and 'lineage'
# duration <- duration %>%
#   separate(
#     col = country_label,
#     into = c("country", "lineage"),
#     sep = "-"
#   )
# 
# # Clean the 'lineage' column by removing text in parentheses
# duration$lineage <- gsub("\\s*\\(.*\\)", "", duration$lineage)
# 
# # Convert 'start_date' to Date format
# duration$start_date <- as.Date(duration$start_date)
# 
# # Calculate 'end_date' based on 'start_date' and 'duration'
# duration$end_date <- duration$start_date + duration$duration
# 
# # --- Step 2: Generate Summary Table ---
# summary_table <- duration %>%
#   group_by(country) %>%
#   summarise(
#     earliest_start_date = min(start_date, na.rm = TRUE),
#     latest_end_date = max(end_date, na.rm = TRUE),
#     unique_lineages = paste(unique(lineage), collapse = ", ")
#   )
# 
# # --- Step 3: Export to Excel ---
# output_file <- "country_summary.xlsx"
# write.xlsx(
#   summary_table,
#   file = output_file,
#   sheetName = "Summary",
#   rowNames = FALSE
# )


# Clean up objects starting with "Nov15" or "Nov25"
rm(list = ls(pattern = "^file"))




country_based_OTC_model <- function(country_name, 
                                    list, 
                                    select = "days_30", 
                                    dominance, 
                                    duration) {
  
  suppressWarnings({
    
    # country_name = "South Korea"; list = gr.list; select = "days_30"; # for test
    # dominance = dominance; duration = duration # for test
    
    gr.days <- list[[select]]
    if (is.null(gr.days)) stop("Error: Selected element does not exist in the list.")
    gr.days <- gr.days[grepl(country_name, names(gr.days))]
    
    if (length(gr.days) == 0) stop("Error: No matching country data found in gr.days.")
    
    # Convert to the OTCs
    sharing.gr <- gr.days
    
    for (i in names(gr.days)) {
      data <- gr.days[[i]]
      if (is.null(data$smoothed_sharing_growth)) next
      sharing.gr[[i]] <- c(i, data$smoothed_sharing_growth[-1])
    }
    
    process_gr_data <- function(gr_list) {
      gr_df <- do.call(rbind, lapply(gr_list, function(x) t(data.frame(x))))
      gr_df <- as.data.frame(gr_df)
      gr_df <- gr_df %>% mutate(across(-1, as.numeric))
      colnames(gr_df)[1] <- "country_label"
      return(gr_df)
    }
    
    sharing_gr_df <- process_gr_data(sharing.gr)
    
    grid.span <- seq(from = 0, to = .4, length.out = 201)
    sharing_otc <- matrix(nrow = nrow(sharing_gr_df), ncol = length(grid.span))
    
    compute_otc <- function(series, grid) {
      sapply(grid, function(x) sum(series > x) / length(series>0))
    }
    
    for (i in seq_len(nrow(sharing_gr_df))) {
      series <- as.numeric(sharing_gr_df[i, -1])
      if (any(is.na(series))) next
      sharing_otc[i, ] <- compute_otc(series, grid.span)
    }
    
    if (!"country_label" %in% colnames(dominance)) {
      stop("Error: 'country_label' column not found in 'dominance'")
    }
    
    dominance_otc <- dominance %>% dplyr::filter(country_label %in% sharing_gr_df$country_label)
    duration_otc <- duration %>% dplyr::filter(country_label %in% sharing_gr_df$country_label)
    
    fused_lasso <- fusedlasso1d(
      y = dominance_otc$mean_share,    
      X = sharing_otc
    )
    
    best_lambda_index <- which.min(fused_lasso$lambda)
    best_lambda <- fused_lasso$lambda[best_lambda_index]
    
    plot_coefficients <- coef(fused_lasso, lambda = .8)$beta
    plot(
      x = grid.span, 
      y = plot_coefficients, 
      type = "b",
      pch = 16,
      cex = 0.5,
      main = bquote(bold(.(country_name) ~ "(" * .(select) * ")")),
      xlab = "Growth Rate",
      ylab = "Coefficients (β)"
    )
    
    plot_result <- recordPlot()
    
    best_coefficients <- coef(fused_lasso, lambda = .05)$beta
    predictions <- sharing_otc %*% best_coefficients
    
    Rsq <- sqrt(cor(dominance_otc$mean_share, predictions)) %>% as.vector()
    
  }) # End suppressWarnings
  
  return(list(
    plot = plot_result,
    plot_ceof = plot_coefficients, 
    prediction_Rsq = Rsq
  ))
}


run_country_OTC_models <- function(country_name) {

    if (!exists("gr.list") || !exists("dominance") || !exists("duration")) {
    stop("Error: 'gr.list', 'dominance', or 'duration' is not defined in the global environment.")
  }
  
  selection_seq <- names(gr.list)
  
  Rsq_seq <- numeric(length(selection_seq))
  plot_list <- vector("list", length(selection_seq))
  plot_coef_list <- vector("list", length(selection_seq))
  
  for (i in seq_along(selection_seq)) {
    current_select <- selection_seq[i]
    cat("Running model for select =", current_select, "\n")
    
    model_result <- country_based_OTC_model(
      country_name = country_name,
      list = gr.list,
      select = current_select,
      dominance = dominance,
      duration = duration
    )
    
    Rsq_seq[i] <- model_result$prediction_Rsq
    plot_list[[i]] <- model_result$plot
    plot_coef_list[[i]] <- model_result$plot_ceof
  }
  
  names(Rsq_seq) <- selection_seq
  names(plot_list) <- selection_seq
  
  min_index <- which.min(Rsq_seq)
  min_Rsq <- Rsq_seq[min_index]
  min_plot <- plot_list[[min_index]]
  min_plot_coef <- plot_coef_list[[min_index]]
  best_select <- selection_seq[min_index]
  
  cat("Minimum prediction_Rsq:", min_Rsq, "\n")
  cat("Best select value:", best_select, "\n")
  
  result <- list(
    Rsq_sequence = Rsq_seq,
    plots = plot_list,
    best_plot = min_plot,
    best_coef = min_plot_coef, 
    best_select = best_select,
    min_Rsq = min_Rsq
  )
  
  return(result)
}


South_Korea    <- run_country_OTC_models(country_name = "South Korea")
United_Kingdom <- run_country_OTC_models(country_name = "United Kingdom")
USA            <- run_country_OTC_models(country_name = "USA")
Australia      <- run_country_OTC_models(country_name = "Australia")
Brazil         <- run_country_OTC_models(country_name = "Brazil")
Canada         <- run_country_OTC_models(country_name = "Canada")


pdf("Manuscripts/plots/Rsq_Sequence.pdf", width = 10, height = 8)

par(mfrow = c(3, 2), mar = c(4, 4, 2, 1), oma = c(2, 2, 4, 2))

plot_Rsq_sequence <- function(country_data, country_name) {
  
  # country_data <- South_Korea; country_name = "South Korea"
  
  highest_index <- which.max(country_data$Rsq_sequence) + 29
  
  plot(
    x = 30 - 1 + seq_along(country_data$Rsq_sequence),         
    y = country_data$Rsq_sequence,                    
    type = "b",                     
    pch = 1,                        
    lty = 1,                        
    lwd = 1.5,                      
    cex = 1,                        
    main = paste0(country_name, "(", highest_index, ")"),  
    xlab = "Days to Input",      
    ylab = "Prediction Rsq"
  )
  
  max_index <- which.max(country_data$Rsq_sequence)
  max_x <- 30 - 1 + max_index
  max_y <- country_data$Rsq_sequence[max_index]
  
  points(max_x, max_y, col = "red", pch = 19, cex = 1.5)
  text(max_x, max_y, labels = round(max_y, 2), pos = 3, col = "red", cex = 0.8)
}

plot_Rsq_sequence(South_Korea, "South Korea")
plot_Rsq_sequence(United_Kingdom, "United Kingdom")
plot_Rsq_sequence(USA, "United States")
plot_Rsq_sequence(Australia, "Australia")
plot_Rsq_sequence(Brazil, "Brazil")
plot_Rsq_sequence(Canada, "Canada")

mtext(
  "R Squared Sequence of Days to Input",
  outer = TRUE,
  cex = 1.5,
  font = 2,
  line = 1
)

dev.off()

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(5, 4, 4, 2) + 0.1)








library(grid)
library(gridExtra)


grid.span <- seq(from = 0, to = 0.4, length.out = 201)

pdf("Manuscripts/plots/Coefficients.pdf", width = 10, height = 8)

par(mfrow = c(3, 2), mar = c(4, 4, 2, 1), oma = c(2, 2, 4, 2))

plot_coef <- function(country_data, country_name) {
  plot(
    x = grid.span, 
    y = country_data$best_coef+.01, 
    type = "b",
    pch = 16,
    cex = 0.5,
    main = bquote(bold(.(country_name))),
    xlab = "Growth Rate",
    ylab = "Coefficients (β)"
  )
}

plot_coef(South_Korea, "South Korea")
plot_coef(United_Kingdom, "United Kingdom")
plot_coef(USA, "United States")
plot_coef(Australia, "Australia")
plot_coef(Brazil, "Brazil")
plot_coef(Canada, "Canada")


mtext(
  "Coefficients for OTCs",
  outer = TRUE,
  cex = 1.5,
  font = 2,
  line = 1
)

dev.off()

par(mfrow = c(1, 1), oma = c(0, 0, 0, 0), mar = c(5, 4, 4, 2) + 0.1)




