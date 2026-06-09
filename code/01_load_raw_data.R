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
source("code/funtions_data_io.R")

# load the raw data
url <- "reference/UKHSA-UConn-variant-modelling/variant_modelling/data/summary_GISAID_20240918.csv"
data <- read.csv(url)

# Extract top 15 countries by total infections
main_countries <- data %>%
  group_by(country) %>%
  summarise(num_total = sum(numerator)) %>%
  arrange(desc(num_total)) %>%
  head(15) %>%
  pull(country)

# Filter data for main countries
data2 <- data %>%
  dplyr::filter(country %in% main_countries) %>% 
  mutate(country = dplyr::recode(country, "USA" = "United States"))


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

target_order <- c("United States", "United Kingdom", "Germany", "Japan", "Denmark",
                  "Canada", "France", "Sweden", "Spain", "Brazil",
                  "Australia", "South Korea", "India", "Italy", "Austria")

country_list <- factor(country_list, levels = target_order)
country_list <- country_list[order(as.numeric(country_list))]

# plot to show the trimmed period
par(mfrow = c(3, 5), mar = c(4, 4, 2, 1), oma = c(2, 2, 3, 1))

for (cty in country_list) {
  df_cty <- data2[data2$country == cty, ] %>%
    dplyr::select(country, date, denominator) %>%
    dplyr::distinct() %>%
    dplyr::arrange(date)

  df_cty$date <- as.Date(df_cty$date)

  df_trim <- trimmed_data_list[[cty]]
  has_trim <- !is.null(df_trim) && nrow(df_trim) > 0

  plot(df_cty$date, df_cty$denominator,
       type = "p", pch = 16, cex = 0.2,
       xlab = "", ylab = "Denominator",
       main = cty)

  if (has_trim) {
    trim_start <- min(df_trim$date, na.rm = TRUE)
    trim_end <- max(df_trim$date, na.rm = TRUE)

    abline(v = trim_start, col = "blue", lty = 1, lwd = 1)
    abline(v = trim_end, col = "blue", lty = 1, lwd = 1)
  }
}

## Sep 25, 2025
# par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
#
# for (cty in c("United States", "South Korea")) {
#   df_cty <- data2[data2$country == cty, ] %>%
#     dplyr::select(country, date, denominator) %>%
#     dplyr::distinct() %>%
#     dplyr::arrange(date)
#
#   df_cty$date <- as.Date(df_cty$date)
#   df_trim <- trimmed_data_list[[cty]]
#   has_trim <- !is.null(df_trim) && nrow(df_trim) > 0
#
#   plot(df_cty$date, df_cty$denominator,
#        type = "n", xlab = "", ylab = "Total COVID-19 cases",
#        main = cty, xaxt = "n", xlim = as.Date(c("2020-01-01", "2024-09-30")))
#
#   segments(df_cty$date, 0, df_cty$date, df_cty$denominator, col = "black", lwd = 0.5)
#
#   axis.Date(1, at = seq(as.Date("2021-01-01"), as.Date("2024-09-30"), by = "3 months"),
#             format = "%b %Y")
#
#   if (has_trim) {
#     trim_start <- min(df_trim$date, na.rm = TRUE)
#     trim_end <- max(df_trim$date, na.rm = TRUE)
#
#     col_trim <- ifelse(cty == "United States", "blue", "orange")
#
#     abline(v = trim_start, col = col_trim, lty = 1, lwd = 1)
#     abline(v = trim_end, col = col_trim, lty = 1, lwd = 1)
#
#     text(trim_start, max(df_cty$denominator, na.rm = TRUE),
#          labels = format(trim_start, "%b %Y"), pos = 4, col = col_trim, cex = 0.7, offset = 0.5)
#     text(trim_end, max(df_cty$denominator, na.rm = TRUE),
#          labels = format(trim_end, "%b %Y"), pos = 2, col = col_trim, cex = 0.7, offset = 0.5)
#   }
# }

data2 <- data2 %>% 
  mutate(share = numerator/denominator,
         date = as.Date(date))

# run once only
# res <- build_lineage_periods(data2, WIN = 14L, SUM_MIN = 10L, NZ_MIN = 3L)
# saveRDS(res, "code/variant_analytical_periods.rds")

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

invalid_combos <- lineage_period_df %>%
  filter(!fully_within_analytical) %>%
  select(country, lineage) %>%
  distinct()

data3 <- data2 %>%
  semi_join(valid_combos, by = c("country", "lineage"))

data_dropped <- data2 %>%
  semi_join(invalid_combos, by = c("country", "lineage"))

data3 <- data3 %>%
  left_join(analytical_period_df, by = "country") %>%
  filter(date >= analytical_start & date <= analytical_end) %>%
  select(-analytical_start, -analytical_end) %>% 
  arrange(country, lineage, date) 

# out_list <- list(
#   data3        = data3,
#   data_dropped = data_dropped
# )
# 
# saveRDS(out_list, file = "code/data3_and_dropped.rds")

data_dropped$country %>% unique()

#########################################################
# to show the dropped variants 
plot_top_lineages_with_periods(
  country_name = "South Korea",
  df = data_dropped,
  analytical_period_df = analytical_period_df,
  n = 6,
  start_date = as.Date("2020-01-01"),
  end_date   = as.Date("2024-08-31")
)

# df_check <- data_dropped %>% 
#   dplyr::filter(country == "United States") %>% 
#   group_by(date) %>% 
#   summarise(sum_share = sum(share), .groups = "drop")
# 
# plot(df_check, type = "l")
