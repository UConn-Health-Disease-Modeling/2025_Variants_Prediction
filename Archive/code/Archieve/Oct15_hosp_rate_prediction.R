# ---- Load Required Libraries ----
library(dplyr)
library(tidyr)
library(mgcv)
library(openxlsx)
library(purrr)
library(fda)
library(forecast)
library(vars)
library(fuzzyjoin)
library(stringr)

################################################################################
# ---- Process Infection-Hospitalization Data ----

# File paths
hosp_url <- "Data/weekly-hospital-admissions-covid-per-million.csv"
daily_est_infec_url <- "Data/daily-new-estimated-infections-of-covid-19.xlsx"

# Load data
hospitalization_data <- read.csv(hosp_url)
daily_estimated_infections <- openxlsx::read.xlsx(daily_est_infec_url)

# Rename columns
colnames(hospitalization_data) <- c("country", "Code", "date", "hosp_per_million")
colnames(daily_estimated_infections) <- c(
  "country", "Code", "date", "ICL_mean", "IHME_mean", "LSHTM_median", "YG_mean", "avg_7_days"
)

# Convert dates
daily_estimated_infections$date <- as.Date(daily_estimated_infections$date, origin = "1899-12-31")
hospitalization_data$date <- as.Date(hospitalization_data$date)

main_countries <- readRDS("Code/ProcessedData/main_countries.rds")

# Combine data
infection_hosp_data <- hospitalization_data %>%
  dplyr::left_join(daily_estimated_infections, by = c("country", "date")) %>%
  dplyr::select(country, date, hosp_per_million, avg_7_days) %>%
  dplyr::mutate(country = ifelse(country == "United States", "USA", country)) %>% 
  dplyr::filter(country %in% main_countries) %>%
  dplyr::mutate(hosp_infection_rate = hosp_per_million / avg_7_days)

infection_hosp_data <- infection_hosp_data %>% dplyr::filter(hosp_infection_rate != Inf)



################################################################################
# ---- Process Variants Data ----

# File paths
variant_url <- "UKHSA-UConn-variant-modelling/variant_modelling/data/summary_GISAID_20240918.csv"
classification_save_url <- "Code/ProcessedData/global_classification_to_append.rds"

# Load data
variant_data <- read.csv(variant_url)
classification_to_append <- readRDS(classification_save_url)

# Clean and transform variant data
variant_data_cleaned <- variant_data %>%
  dplyr::filter(country %in% main_countries) %>%
  dplyr::mutate(country = ifelse(country == "USA", "United States", country)) %>%
  dplyr::arrange(country, date) %>%
  dplyr::rename(n = numerator) %>%
  dplyr::select(country, date, lineage, n)

# Classify variants
classified_data <- variant_data_cleaned %>%
  dplyr::mutate(decimal_lineage = paste0(lineage, ".")) %>%
  fuzzyjoin::fuzzy_left_join(
    classification_to_append,
    by = "decimal_lineage",
    match_fun = stringr::str_starts
  ) %>%
  dplyr::filter(
    is.na(classified_unasliased) |
      stringr::str_length(classified_unasliased) == max(stringr::str_length(classified_unasliased)),
    .by = c(lineage, date)
  ) %>%
  dplyr::mutate(
    length_class = stringr::str_count(classified_unasliased, "\\."),
    length_lineage = stringr::str_count(lineage, "\\."),
    classified_label = dplyr::if_else(
      (length_lineage - length_class) <= 1,
      classified_label,
      NA_character_
    ),
    classified_label = dplyr::if_else(
      is.na(classified_label),
      "Other",
      classified_label
    )
  ) %>%
  dplyr::filter(classified_label != "Other") %>%
  dplyr::select(country, date, lineage, n, classified_label) %>%
  dplyr::group_by(country, classified_label, date) %>%
  dplyr::summarise(n_total = sum(n, na.rm = TRUE), .groups = "drop")

# ---- Complete Dates to be Continuous ----
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
    tidyr::replace_na(list(n_total = 0)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-min_date, -max_date)  # Drop temporary columns
}

classified_data$date <- as.Date(classified_data$date)

classified_data <- classified_data %>%
  complete_dates() %>%
  dplyr::rename(numerator = n_total) %>%
  dplyr::group_by(country, date) %>%
  dplyr::mutate(denominator = sum(numerator, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(share = numerator / denominator)

classified_data$country[classified_data$country == "United States"] <- "USA"

# ---- Process Variants for Each Country ----
process_variants <- function(data, country_name) {
  variants <- data %>% dplyr::filter(country == country_name)
  dates <- unique(variants$date)
  variants_list <- lapply(unique(variants$classified_label), function(label) {
    variants %>% dplyr::filter(classified_label == label) %>% dplyr::pull(share)
  })
  names(variants_list) <- unique(variants$classified_label)
  variants_df <- cbind(date = dates, as.data.frame(variants_list))
  return(variants_df)
}

# ---- Process Variants for Each Country ----

variants_share_list <- list()
unique_countries <- unique(infection_hosp_data$country)

for (country_name in unique_countries) {
  variants_share_list[[country_name]] <- process_variants(classified_data, country_name)
}

combined_data_list <- list()
unique_countries <- unique(infection_hosp_data$country)

for (country_name in unique_countries) {
  
  # country_name = "United Kingdom"
  
  hosp_data <- infection_hosp_data %>%
    dplyr::filter(country == country_name) %>%
    dplyr::select(date, hosp_infection_rate)
  
  combined_data <- variants_share_list[[country_name]] %>%
    dplyr::left_join(hosp_data, by = "date") %>%
    dplyr::filter(!is.na(hosp_infection_rate))
  
  combined_data <- combined_data[, apply(combined_data, 2, function(col) !all(col == 0))]
  
  combined_data_list[[country_name]] <- combined_data
  
}




library(forecast)
library(Metrics)

pdf("Manuscripts/plots/Hosp_Infection_Rate_Plots.pdf", width = 8, height = 6)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

coefficients_list <- list()

Top_4.eg <- c("South Korea", "United Kingdom", "USA", "Spain")

for (country in Top_4.eg) {
  
  # country = "United Kingdom"
  
  data.eg <- combined_data_list[[country]]
  
  y <- data.eg$hosp_infection_rate
  predictors <- as.matrix(data.eg[, -c(1, ncol(data.eg))])
  
  arimax_model <- auto.arima(y, xreg = predictors)
  fitted_values <- fitted(arimax_model)
  mape_value <- mape(y, fitted_values)
  
  cat("Country:", country, " - MAPE:", round(mape_value * 100, 2), "%\n")
  
  # Save coefficients as a dataframe
  coefficients <- coef(arimax_model)
  coefficients_df <- data.frame(
    Name = names(coefficients),
    Coefficient = as.numeric(coefficients)
  )
  coefficients_list[[country]] <- coefficients_df
  
  # Plot actual vs fitted values
  plot(data.eg$date, y, type = "l", col = "blue", lwd = 2,
       main = country, xlab = "", ylab = "Hosp-Infec Rate", xaxt = "n")
  lines(data.eg$date, fitted_values, col = "red", lwd = 2)
  axis.Date(1, at = seq(min(data.eg$date), max(data.eg$date), by = "6 months"), format = "%b-%Y")
  legend("topright", legend = c("Actual", "Fitted"), col = c("blue", "red"),
         lty = 1, lwd = 2, cex = 0.8, bty = "n")
}

par(mfrow = c(1, 1))
dev.off()


file_paths <- list(
  dominance = "Code/ProcessedData/Nov15_variants_dominance.rds",
  duration  = "Code/ProcessedData/Nov25_variants_duration.rds"
)

# Load data
dominance <- readRDS(file_paths$dominance)

dominance <- dominance %>%
  separate(col = country_label, into = c("country", "lineage"), sep = "-") %>% 
  dplyr::filter(country %in% c("South Korea", "United Kingdom", "USA", "Spain"))


coefficients_list$`South Korea`$country <- "South Korea"
coefficients_list$`United Kingdom`$country <- "United Kingdom"
coefficients_list$USA$country <- "USA"
coefficients_list$Spain$country <- "Spain"


coef_df <- rbind(coefficients_list$`South Korea`, 
                 coefficients_list$`United Kingdom`, 
                 coefficients_list$USA, 
                 coefficients_list$Spain)

colnames(coef_df)[1] <- "lineage"

coef_lineage <- unique(coef_df$lineage)
dominance_lineage <- unique(dominance$lineage)


mapping <- data.frame(
  original = coef_lineage,
  matched = sapply(coef_lineage, function(x) {
    dominance_lineage[which.min(stringdist::stringdist(x, dominance_lineage, method = "jw"))]
  })
)

mapping$matched[11] <- "BA.1 (BA.1)"
mapping$matched[14] <- "BA.2 (BA.2)"
mapping$matched[16] <- "BA.5 (BA.5)"


coef_df <- coef_df %>%
  left_join(mapping, by = c("lineage" = "original")) %>%
  mutate(lineage = matched) %>%
  dplyr::select(-matched)

comb_df <- coef_df %>% left_join(dominance, by = c("country", "lineage")) %>% 
  dplyr::filter(!is.na(mean_share))

comb_df$lineage[comb_df$lineage == "BA.1 (BA.1)"] <- "BA.1"
comb_df$lineage[comb_df$lineage == "BA.1.1 (BA.1.1)"] <- "BA.1.1"
comb_df$lineage[comb_df$lineage == "BA.1.17 (BA.1.17)"] <- "BA.1.17"
comb_df$lineage[comb_df$lineage == "BA.2 (BA.2)"] <- "BA.2"
comb_df$lineage[comb_df$lineage == "BA.4 (BA.4)"] <- "BA.4"
comb_df$lineage[comb_df$lineage == "BA.5 (BA.5)"] <- "BA.5"
comb_df$lineage[comb_df$lineage == "BA.5.2 (BA.5.2)"] <- "BA.5.2"
comb_df$lineage[comb_df$lineage == "BA.1.17 (BA.1.17)"] <- "BA.1.17"

comb_df.uk <- comb_df %>% dplyr::filter(country == "United Kingdom")

y_range <- range(comb_df.uk$Coefficient)
y_range <- c(y_range[1] - 0.01, y_range[2] + 0.01) 

plot(comb_df.uk$mean_share, comb_df.uk$Coefficient, 
     xlab = "Dominance", 
     ylab = "Hosp Coefficient",
     main = "United Kingdom",
     pch = 19, col = "blue",
     ylim = y_range, 
     cex.main = 1.2, font.main = 2)  


text(comb_df.uk$mean_share, comb_df.uk$Coefficient, 
     labels = comb_df.uk$lineage, pos = 1, cex = 0.7, col = "black")





comb_df.us <- comb_df %>% dplyr::filter(country == "USA")

y_range <- range(comb_df.us$Coefficient)
y_range <- c(y_range[1] - 0.01, y_range[2] + 0.01) 

plot(comb_df.us$mean_share, comb_df.us$Coefficient, 
     xlab = "Dominance", 
     ylab = "Hosp Coefficient",
     main = "United States",
     pch = 19, col = "blue",
     ylim = y_range, 
     cex.main = 1.2, font.main = 2)  


text(comb_df.us$mean_share, comb_df.us$Coefficient, 
     labels = comb_df.us$lineage, pos = 1, cex = 0.7, col = "black")



comb_df.sk <- comb_df %>% dplyr::filter(country == "South Korea")

y_range <- range(comb_df.sk$Coefficient)
y_range <- c(y_range[1] - 0.01, y_range[2] + 0.01) 

plot(comb_df.sk$mean_share, comb_df.sk$Coefficient, 
     xlab = "Dominance", 
     ylab = "Hosp Coefficient",
     main = "South Korea",
     pch = 19, col = "blue",
     ylim = y_range, 
     cex.main = 1.2, font.main = 2)  


text(comb_df.sk$mean_share, comb_df.sk$Coefficient, 
     labels = comb_df.sk$lineage, pos = 1, cex = 0.7, col = "black")


