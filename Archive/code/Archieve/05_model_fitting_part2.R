# ------------------------------------------------------------------------------
rm(list = ls())

# Load necessary packages
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(lubridate)
library(stringr)

source("code/functions_ml_utils.R")

data <- readRDS( "code/data3_and_dropped.rds")[[1]]

if (!inherits(data$date, "Date")) data$date <- as.Date(data$date)

thr <- 0.01  # 1% threshold

summary_tbl <- data %>%
  dplyr::group_by(country, lineage) %>%
  dplyr::group_modify(~{
    df <- dplyr::arrange(.x, date)
    s  <- df$share
    d  <- df$date
    
    idx_reach <- which(is.finite(s) & s >= thr)
    
    if (length(idx_reach) == 0) {
      tibble::tibble(
        first_reach_date      = as.Date(NA),
        status                = "excluded",
        reason                = "never_reached_1pct",
        duration_days_ge_1pct = NA_integer_,
        duration_obs_ge_1pct  = 0L
      )
    } else {
      i0 <- idx_reach[1]                                   # first reach
      seg <- is.finite(s[i0:length(s)]) & s[i0:length(s)] >= thr
      first_break <- which(!seg)[1]
      consec_len  <- if (is.na(first_break)) length(seg) else (first_break - 1L)
      
      i1 <- i0 + consec_len - 1L                           # end index
      start_date <- d[i0]
      end_date   <- d[i1]
      dur_days   <- as.integer(end_date - start_date) + 1L
      
      tibble::tibble(
        first_reach_date      = start_date,
        status                = "included",
        reason                = NA_character_,
        duration_days_ge_1pct = dur_days,                   # calendar days
        duration_obs_ge_1pct  = consec_len                  # consecutive observations
      )
    }
  }) %>%
  dplyr::ungroup() 

summary_tbl <- summary_tbl %>% 
  dplyr::mutate(duration_days_ge_1pct = ifelse(is.na(duration_days_ge_1pct), 0, duration_days_ge_1pct))


eligible_data <- summary_tbl %>% dplyr::filter(duration_days_ge_1pct > 35)

eligible_data %>% filter(country == "United States")
table(eligible_data$country)




