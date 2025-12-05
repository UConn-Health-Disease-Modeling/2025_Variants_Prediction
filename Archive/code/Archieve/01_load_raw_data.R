# clear the environment
rm(list = ls())

# Load necessary packages
library(dplyr)
library(purrr)
library(lubridate)
library(stringr)
library(zoo)

source("reference/reference code/growth_rate_function_original.R")
source("reference/reference code/functions_0706.R")

# ------------------------------------------------------------------------------
# Processing Stages: 
# 1) Select Top 15 main countries and get the analytical period for each country 
# 2) Append classified alias 
# 3) Re-define the epidemic period of the variant. 
# 4) Define the measurements: "smoothed maximum" and "duration"

# Justification of 1):
# To ensure the reliability of variant share estimates, we first identified 
# analytically valid periods for each country based on national-level sequencing 
# coverage following World Health Organization (WHO) guidelines. Specifically, 
# a valid analytical period was defined as a continuous stretch of at least three 
# calendar weeks in which each week contained a minimum of 100 sequences and at 
# least four days with non-zero sampling. This standard reflects WHO recommendations 
# that reliable variant prevalence estimation requires consistent and representative 
# sequencing, with at least 100 samples per week per geographic area (WHO, 2021). 
# It also aligns with empirical findings from Chen et al. (2022), who showed that 
# sparse or irregular early sampling leads to considerable noise and bias in 
# share-based modeling. We then retained only those variant–country combinations 
# whose entire observed share trajectories fell completely within these valid 
# national-level periods. This filtering strategy ensures that both early growth 
# patterns and peak share estimates are derived from statistically meaningful data 
# and are not distorted by low-coverage or episodic sampling.

# Justification of 3): 
# Du Plessis, L., McCrone, J. T., Zarebski, A. E., Hill, V., Ruis, C., Gutierrez, B., … & Pybus, O. G. (2021). 
# Establishment and lineage dynamics of the SARS-CoV-2 epidemic in the UK. Science, 371(6530), 708–712. 
# https://doi.org/10.1126/science.abf2946
# "Lineages were considered established in the UK if they resulted in at least 10 sequenced cases over a period of 14 days or longer…"
# 
# Justification of 4): 
# Volz, E., Mishra, S., Chand, M., Barrett, J. C., Johnson, R., Geidelberg, L., et al. (2021).
# Assessing transmissibility of SARS-CoV-2 lineage B.1.1.7 in England.
# Nature, 593(7858), 266–269.
# https://doi.org/10.1038/s41586-021-03470-x
# "Used smoothed variant share curves to estimate variant dynamics and growth."
# 
# Althaus, C. L., et al. (2021). Rapid spread of the SARS-CoV-2 variant of concern 202012/01 in Switzerland.
# medRxiv. https://doi.org/10.1101/2021.01.10.21249380
# "Applied a LOESS smoother to daily frequency estimates."
# 
# UK Health Security Agency (UKHSA) Technical Briefings
# Frequently use 7-day rolling averages and smoothed curves to report peak share and variant dynamics.
# 
# “We defined the peak share of a variant as the maximum of a 7-day moving average of its daily share, 
# consistent with methods used in previous genomic surveillance studies (e.g., Volz et al., 2021; Althaus et al., 2021). 
# This approach reduces noise from day-to-day fluctuations and sequencing variation.”


# load the raw data
url <- "reference/UKHSA-UConn-variant-modelling/variant_modelling/data/summary_GISAID_20240918.csv"
data <- read.csv(url)

# ranked_country <- data %>% 
#   dplyr::group_by(country) %>% 
#   dplyr::summarise(total_case = sum(numerator)) %>% 
#   dplyr::arrange(desc(total_case)) %>% 
#   dplyr::slice_head(n = 15)
# 
# # list.files()
# saveRDS(ranked_country, "code/ranked_country.rds")

# Extract top 15 countries by total infections
main_countries <- data %>%
  group_by(country) %>%
  summarise(num_total = sum(numerator)) %>%
  arrange(desc(num_total)) %>%
  head(15) %>%
  pull(country)

# Filter data for main countries
data2 <- data %>%
  dplyr::filter(country %in% main_countries)


country_list <- unique(data2$country)

trimmed_data_list <- list()
for (cty in country_list) {
  df_cty <- data2[data2$country == cty, ] %>%
    dplyr::select(country, date, denominator) %>%
    dplyr::distinct() %>%
    dplyr::arrange(date)
  
  df_cty$date <- as.Date(df_cty$date)
  
  df_cty <- df_cty %>%
    mutate(week = floor_date(date, unit = "week", week_start = 1))
  
  weekly_summary <- df_cty %>%
    group_by(week) %>%
    summarise(
      n_total = sum(denominator, na.rm = TRUE),
      n_days_ge20 = sum(denominator >= 20, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(week_valid = n_total >= 100 & n_days_ge20 >= 5)
  
  weekly_summary <- weekly_summary %>%
    mutate(run_id = data.table::rleid(week_valid))
  
  valid_runs <- weekly_summary %>%
    group_by(run_id) %>%
    summarise(
      start_week = min(week),
      end_week = max(week),
      run_length = n(),
      valid = first(week_valid),
      .groups = "drop"
    ) %>%
    filter(valid == TRUE, run_length >= 3)
  
  if (nrow(valid_runs) > 0) {
    overall_start <- min(valid_runs$start_week)
    overall_end <- max(valid_runs$end_week) + days(6)
    
    trimmed_df <- df_cty %>%
      filter(date >= overall_start & date <= overall_end)
    
    trimmed_data_list[[cty]] <- trimmed_df
  } else {
    message(paste("No WHO-valid period found for", cty))
  }
}

# # plot to show the trimmed period
# par(mfrow = c(3, 5), mar = c(4, 4, 2, 1), oma = c(2, 2, 3, 1))
# 
# for (cty in country_list) {
#   df_cty <- data2[data2$country == cty, ] %>%
#     dplyr::select(country, date, denominator) %>%
#     dplyr::distinct() %>%
#     dplyr::arrange(date)
# 
#   df_cty$date <- as.Date(df_cty$date)
# 
#   df_trim <- trimmed_data_list[[cty]]
#   has_trim <- !is.null(df_trim) && nrow(df_trim) > 0
# 
#   plot(df_cty$date, df_cty$denominator,
#        type = "p", pch = 16, cex = 0.2,
#        xlab = "", ylab = "Denominator",
#        main = cty)
# 
#   # 如果该国家有裁剪区间，则添加蓝色竖线
#   if (has_trim) {
#     trim_start <- min(df_trim$date, na.rm = TRUE)
#     trim_end <- max(df_trim$date, na.rm = TRUE)
# 
#     abline(v = trim_start, col = "blue", lty = 1, lwd = 1)
#     abline(v = trim_end, col = "blue", lty = 1, lwd = 1)
#   }
# }

data2 <- data2 %>% 
  mutate(share = numerator/denominator,
         date = as.Date(date))

lineage_periods <- list()
combo_list <- data2 %>%
  distinct(country, lineage)
n_total <- nrow(combo_list)

# # Only run once
# pb <- txtProgressBar(min = 0, max = n_total, style = 3)
# 
# for (i in seq_len(n_total)) {
#   cty <- combo_list$country[i]
#   lin <- combo_list$lineage[i]
#   
#   # 更新进度条
#   setTxtProgressBar(pb, i)
#   
#   df_sub <- data2 %>%
#     filter(country == cty, lineage == lin) %>%
#     arrange(date)
#   
#   if (nrow(df_sub) < 14) next
#   
#   df_sub <- df_sub %>%
#     mutate(
#       roll_sum = zoo::rollapply(numerator, width = 14, FUN = sum, fill = NA, align = "left"),
#       roll_nonzero_days = zoo::rollapply(numerator > 0, width = 14, FUN = sum, fill = NA, align = "left")
#     )
#   
#   valid_rows <- df_sub %>%
#     filter(roll_sum >= 10, roll_nonzero_days >= 3)
#   
#   if (nrow(valid_rows) == 0) next
#   
#   start_date <- min(valid_rows$date)
#   end_date <- max(valid_rows$date) + 13
#   
#   lineage_periods[[paste(cty, lin, sep = "|")]] <- data.frame(
#     country = cty,
#     lineage = lin,
#     start_date = start_date,
#     end_date = end_date
#   )
# }
# 
# close(pb)
# 
# lineage_period_df <- bind_rows(lineage_periods)
# saveRDS(lineage_period_df, "code/variant_analytical_periods.rds")

lineage_period_df <- readRDS("code/variant_analytical_periods.rds")

analytical_period_df <- lapply(trimmed_data_list, function(df) {
  data.frame(
    country = unique(df$country),
    analytical_start = min(df$date),
    analytical_end = max(df$date)
  )
}) %>%
  bind_rows()

lineage_period_df <- lineage_period_df %>%
  left_join(analytical_period_df, by = "country") %>%
  mutate(
    fully_within_analytical = start_date >= analytical_start & end_date <= analytical_end
  )

valid_combos <- lineage_period_df %>%
  filter(fully_within_analytical) %>%
  select(country, lineage) %>%
  distinct()

data3 <- data2 %>%
  semi_join(valid_combos, by = c("country", "lineage"))

data3 <- data3 %>%
  left_join(analytical_period_df, by = "country") %>%
  filter(date >= analytical_start & date <= analytical_end) %>%
  select(-analytical_start, -analytical_end)

# saveRDS(data3, "code/variant_trimmed.rds")



# ############################################
# data2 %>% 
#   dplyr::select(country, lineage) %>% 
#   distinct() %>% 
#   dim()
# 
# 
# data3 %>% 
#   dplyr::select(country, lineage) %>% 
#   distinct() %>% 
#   dplyr::select(country) %>% 
#   unlist() %>% 
#   table()
# 
# 
# 
# 
# 
# data2 %>% head(10)
# data3 %>% head(10)
