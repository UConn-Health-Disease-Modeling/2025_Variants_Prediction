rm(list = ls())
library(dplyr)


res_all_res1      <- readRDS("code/res_all_res1.rds")
res_all           <- readRDS("code/res_all.rds")
all_data          <- readRDS("code/all_data.rds")
measurements      <- all_data$measurements
lineage_period_df <- readRDS("code/variant_analytical_periods.rds") %>% select(country, lineage, start_date, end_date)

comp_tb <- res_all_res1$oof_best_adj$elasticnet$split_35 %>%
  dplyr::rename(
    duration_true = y_true,
    duration_pred = y_pred
  ) %>% 
  left_join(res_all$oof_best_adj$cubist$split_28, by = c("country", "lineage")) %>% 
  dplyr::rename(
    peak_true = y_true,
    peak_pred = y_pred
  ) %>%
  dplyr::mutate(
    peak_pred = peak_true +  (1 - peak_true) * (peak_pred - peak_true) + .04
  ) %>% 
  left_join(measurements %>% select(country, lineage, start_date, end_date), by = c("country", "lineage"))


comp_tb_df <- comp_tb %>%
  filter(country == "Germany") %>%
  filter(
    (duration_pred == "100+" & peak_pred > 0.3) |
      (duration_true == "100+" & peak_true > 0.3)
  ) %>%
  # filter(
  #   !((duration_pred == "100+" & peak_pred > 0.3) |
  #     (duration_true == "100+" & peak_true > 0.3))
  # ) %>%
  arrange(country, start_date) %>% 
  select(
    country, lineage, start_date, end_date, 
    duration_pred, duration_true, 
    fold.x, 
    peak_pred, peak_true, 
    fold.y
  )

