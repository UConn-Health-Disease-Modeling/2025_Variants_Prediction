# clear the environment
rm(list = ls())

# Load necessary packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, ggplot2, lubridate, patchwork, purrr, randomcoloR, tidyr, zoo)


source("code2/funtions_data_io.R")
source("code2/classify_lineage.R")

# Load alias list for lineage classification
source("reference/reference code/functions_0706.R")
alias_map <- get_alias()
data3 <- readRDS("code2/data3_and_dropped.rds")[["data3"]]

classified_data <- process_classified_variants(
  data = data3,
  alias_list = alias_map,
  p_lim = 0.01
)

classified_data <- classified_data %>%
  dplyr::filter(!grepl("Unassigned", country_label)) %>%
  tidyr::separate(country_label, into = c("country", "lineage"), sep = "-")





classified_data <- readRDS("code2/all_data.rds")[["classified_data"]]
# Count the number of variants for each country
classified_data %>%
  dplyr::distinct(country, lineage) %>%
  dplyr::select(country) %>% unlist() %>% table()
#     Australia        Austria         Brazil         Canada        Denmark         France        Germany
#            43             23             41            89              57             72             95
#         India          Italy          Japan    South Korea          Spain         Sweden United Kingdom
#            30             35             83             33             47             47             75
# United States
#           104



# ##############################################################################
# define the measurements

library(dplyr)
library(zoo)

peak_by_variant <- classified_data %>%
  group_by(country, lineage) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    mean_14d = rollapply(
      share, width = 14, FUN = mean,
      align = "right", partial = FALSE, fill = NA, na.rm = TRUE
    )
  ) %>%
  summarise(
    peak_share    = max(mean_14d, na.rm = TRUE),
    peak_date     = date[which.max(mean_14d)],
    second_1pct_date = {
      over_idx <- which(share >= 0.01)
      if (length(over_idx) >= 3) date[over_idx[3]] else as.Date(NA)
    },
    .groups = "drop"
  ) %>%
  filter(!is.na(second_1pct_date)) %>%
  mutate(
    time_to_peak = as.numeric(peak_date - second_1pct_date),
    eps = 0.5 / (nrow(.) + 1),
    peak_share_clamped = pmin(pmax(peak_share, eps), 1 - eps),
    peak_share_logit = qlogis(peak_share_clamped)
  )

# plot_data <- peak_by_variant %>%
#   select(country, lineage, peak_share_3, peak_share_7, peak_share_14, peak_share) %>%
#   pivot_longer(
#     cols = starts_with("peak_share"),
#     names_to = "window",
#     values_to = "value"
#   ) %>%
#   mutate(
#     window = factor(window,
#                     levels = c("peak_share_3", "peak_share_7", "peak_share_14", "peak_share"),
#                     labels = c("3-day", "7-day", "14-day", "30-day")),
#     id = paste(country, lineage, sep = "_")
#   )
# 
# id_order <- plot_data %>%
#   filter(window == "30-day") %>%
#   arrange(desc(value)) %>%
#   pull(id)
# 
# plot_data$id <- factor(plot_data$id, levels = id_order)
# 
# ggplot(plot_data, aes(x = id, y = value, color = window, shape = window)) +
#   geom_point(size = 1.5, alpha = 0.8) +
#   labs(x = NULL, y = "Peak share", color = "Window", shape = "Window") +
#   theme_bw(base_size = 13) +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank()
#   )
# 
# 
# summary_tbl <- tibble(
#   window = c("3-day", "7-day", "14-day"),
#   MAD = c(
#     mean(abs(peak_by_variant$peak_share_3  - peak_by_variant$peak_share), na.rm = TRUE),
#     mean(abs(peak_by_variant$peak_share_7  - peak_by_variant$peak_share), na.rm = TRUE),
#     mean(abs(peak_by_variant$peak_share_14 - peak_by_variant$peak_share), na.rm = TRUE)
#   ),
#   MRD = c(
#     mean(abs(peak_by_variant$peak_share_3  - peak_by_variant$peak_share) / 
#            pmax(peak_by_variant$peak_share, 1e-6), na.rm = TRUE),
#     mean(abs(peak_by_variant$peak_share_7  - peak_by_variant$peak_share) / 
#            pmax(peak_by_variant$peak_share, 1e-6), na.rm = TRUE),
#     mean(abs(peak_by_variant$peak_share_14 - peak_by_variant$peak_share) / 
#            pmax(peak_by_variant$peak_share, 1e-6), na.rm = TRUE)
#   ),
#   Variance_Ratio = c(
#     var(peak_by_variant$peak_share_3,  na.rm = TRUE) / var(peak_by_variant$peak_share, na.rm = TRUE),
#     var(peak_by_variant$peak_share_7,  na.rm = TRUE) / var(peak_by_variant$peak_share, na.rm = TRUE),
#     var(peak_by_variant$peak_share_14, na.rm = TRUE) / var(peak_by_variant$peak_share, na.rm = TRUE)
#   ),
#   Pearson = c(
#     cor(peak_by_variant$peak_share_3,  peak_by_variant$peak_share, use = "complete.obs"),
#     cor(peak_by_variant$peak_share_7,  peak_by_variant$peak_share, use = "complete.obs"),
#     cor(peak_by_variant$peak_share_14, peak_by_variant$peak_share, use = "complete.obs")
#   ),
#   Spearman = c(
#     cor(peak_by_variant$peak_share_3,  peak_by_variant$peak_share, method = "spearman", use = "complete.obs"),
#     cor(peak_by_variant$peak_share_7,  peak_by_variant$peak_share, method = "spearman", use = "complete.obs"),
#     cor(peak_by_variant$peak_share_14, peak_by_variant$peak_share, method = "spearman", use = "complete.obs")
#   ),
#   Kendall = c(
#     cor(peak_by_variant$peak_share_3,  peak_by_variant$peak_share, method = "kendall", use = "complete.obs"),
#     cor(peak_by_variant$peak_share_7,  peak_by_variant$peak_share, method = "kendall", use = "complete.obs"),
#     cor(peak_by_variant$peak_share_14, peak_by_variant$peak_share, method = "kendall", use = "complete.obs")
#   )
# )
#
#
# summary_tbl
# # > summary_tbl
# # # A tibble: 3 × 7
# #   window    MAD   MRD Variance_Ratio Pearson Spearman Kendall
# #   <chr>   <dbl> <dbl>          <dbl>   <dbl>    <dbl>   <dbl>
# # 1 3-day  0.0447 0.446           1.33   0.969    0.982   0.892
# # 2 7-day  0.0230 0.216           1.16   0.989    0.994   0.939
# # 3 14-day 0.0118 0.106           1.08   0.996    0.998   0.966
# 
# # •	Shorter windows increase volatility.
# # The 3-day window shows the largest deviations from 30-day peaks, while 14-day is nearly identical.
# # •	Error decreases with window length.
# # MAD drops from 0.045 → 0.012, and MRD from 45 % → 11 %, showing smoother estimates with longer averaging.
# # •	Variance stabilizes.
# # Variance ratios (1.33 → 1.08) confirm short windows amplify noise, long windows dampen it.
# # •	High consistency overall.
# # Pearson ≥ 0.97, Spearman ≥ 0.98, Kendall ≥ 0.89 indicate strong agreement in both magnitude and rank.
# # •	Conclusion.
# # 7-day and 14-day windows provide nearly the same peak estimates as 30-day, proving the model’s robustness to window choice.






# peak_by_variant <- peak_by_variant %>% 
#   dplyr::select(country, lineage, peak_share, second_1pct_date, peak_date, time_to_peak) %>% 
#   dplyr::filter(time_to_peak >= 35)
# 
# duration_summary <- peak_by_variant %>%
#   group_by(country) %>%
#   summarise(mean_duration_to_peak =mean(time_to_peak, na.rm = TRUE), .groups = "drop")
# duration_summary$mean_duration_to_peak <- round(duration_summary$mean_duration_to_peak)

peak_by_variant$time_to_peak %>% mean() %>% round()

span_stats <- classified_data %>%
  group_by(country, lineage) %>%
  summarise(
    start_date    = min(date, na.rm = TRUE),
    end_date      = max(date, na.rm = TRUE),
    duration_days = as.integer(end_date - start_date) + 1L,
    total_cases   = sum(numerator), 
    days_above_10 = sum(share > 0.05, na.rm = TRUE),
    .groups = "drop"
  )

# # raw distribution is extremely over skewed
# hist(peak_by_variant$peak_share)
# hist(peak_by_variant$peak_share_logit)

span_stats <- span_stats %>%
  mutate(
    days_above_10_cat = case_when(
      days_above_10 == 0 ~ "0",
      days_above_10 <= 30 ~ "1–30",
      days_above_10 <= 100 ~ "31–100",
      TRUE ~ "100+"
    ),
    days_band = factor(days_above_10_cat, levels = c("0","1–30","31–100","100+"))
  )

# # use the categorical labels
# ggplot(span_stats, aes(x = days_band)) +
#   geom_bar(fill = "grey70", color = "grey30") +
#   labs(x = "Days above 10% share (bands)", y = "Count") + theme_bw()

measurements <- span_stats %>%
  left_join(peak_by_variant, by = c("country", "lineage")) %>%
  dplyr::select(country, lineage, 
                start_date, end_date, 
                peak_share, peak_share_logit, 
                days_above_10, days_above_10_cat, 
                duration_days, total_cases, eps, peak_share_clamped)


# measurements$days_above_10_cat %>% table()

# save the plots (TIFF, 600 dpi, LZW) — paths unchanged

first_three <- c("United States","United Kingdom","South Korea")

countries <- classified_data %>%
  dplyr::distinct(country) %>%
  dplyr::pull() %>%
  sort()

legends_list_all <- measurements %>%
  dplyr::filter(days_above_10 > 100, peak_share > 0.3) %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(lineages = list(unique(lineage)), .groups = "drop") %>%
  { setNames(.$lineages, .$country) }

mk_legends <- function(cntries) {
  setNames(lapply(cntries, function(cty) {
    labs <- legends_list_all[[cty]]
    if (is.null(labs) || length(labs) == 0) character(0) else labs
    unique(c(labs, "others"))
  }), cntries)
}

seed <- 6
out_dir <- "result/figs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

res_first <- build_variant_panels(
  data         = classified_data,
  countries    = first_three,
  k            = 5,
  start_date   = as.Date("2020-01-01"),
  end_date     = as.Date("2024-12-31"),
  legends_list = mk_legends(first_three),
  seed         = seed
)

res_first$plot_combined

ggplot2::ggsave(
  filename    = file.path(out_dir, "variant_stack_3x1(1).tiff"),
  plot        = res_first$plot_combined,
  width       = 12,
  height      = 8,
  units       = "in",
  dpi         = 600,
  compression = "lzw"
)

remaining <- setdiff(countries, first_three)
chunks <- split(remaining, ceiling(seq_along(remaining) / 3))

for (i in seq_along(chunks)) {
  trio <- chunks[[i]]
  res_i <- build_variant_panels(
    data         = classified_data,
    countries    = trio,
    k            = 5,
    start_date   = as.Date("2020-01-01"),
    end_date     = as.Date("2024-12-31"),
    legends_list = mk_legends(trio),
    seed         = seed
  )
  
  fn <- file.path(out_dir, sprintf("variant_stack_3x1(%d).tiff", i + 1))
  
  ggplot2::ggsave(
    filename    = fn,
    plot        = res_i$plot_combined,
    width       = 14,
    height      = 9,
    units       = "in",
    dpi         = 600,
    compression = "lzw"
  )
}

# # save the data
# all_data <- list(
#   measurements   = measurements,
#   classified_data = classified_data
# )
#
# saveRDS(all_data, file = "code/all_data.rds")
