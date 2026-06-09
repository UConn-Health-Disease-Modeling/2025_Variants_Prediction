# ------------------------------------------------------------------------------
# Reproducible build of early-period feature sets (day14/day21/day28/day35)
# for the top-15 countries, then save as a single RDS list.
# ------------------------------------------------------------------------------

rm(list = ls())
library(dplyr)
library(stringr)
source("code2/funtions_data_io.R")   # defines prepare_model_data(), extract_features_from_model_data()

all_data       <- readRDS("code2/all_data.rds")

measurements    <- all_data$measurements
classified_data <- all_data$classified_data %>% 
  dplyr::filter(
    lineage %in% measurements$lineage
  )

# # table(measurements$days_above_10_cat)
# 
# # classified_data_us <- classified_data %>% 
# #   filter(country == "United States") %>% 
# #   filter(lineage == "JN.1")
# # 
# # # classified_data_us <- classified_data_us[-(1:935), ]
# # 
# # lin_name <- unique(classified_data_us$lineage)
# # cty_name <- unique(classified_data_us$country)
# # 
# # plot(classified_data_us$date, classified_data_us$share,
# #      type = "p", ylim = c(0, 1),
# #      xlab = "", ylab = "Share (%)",
# #      pch = 18,
# #      cex = 0.5,
# #      yaxt = "n",
# #      xaxt = "n",  
# #      main = paste(cty_name, lin_name, sep = " - "),
# #      cex.main = 1.5)
# # 
# # axis(2, at = seq(0, 1, by = 0.1), labels = paste0(seq(0, 100, by = 10), "%"))
# # 
# # axis.Date(1,
# #           at = seq(as.Date("2021-01-01"), as.Date("2024-09-30"), by = "6 months"),
# #           format = "%b %Y")
# # 
# # idx_max <- which.max(classified_data_us$share)
# # x_max <- classified_data_us$date[idx_max]
# # y_max <- classified_data_us$share[idx_max]
# # 
# # points(x_max, y_max, col = "red", pch = 16, cex = .8)
# # text(x_max, y_max + 0.05, paste0(round(y_max * 100, 1), "%"),
# #      col = "red", cex = 0.9)
# # 
# # above_10 <- classified_data_us %>%
# #   filter(share > 0.10)
# # 
# # if (nrow(above_10) > 0) {
# #   block_start <- min(above_10$date, na.rm = TRUE)
# #   block_end <- max(above_10$date, na.rm = TRUE)
# #   segments(block_start, -0.02, block_end, -0.02, col = "blue", lwd = 5, xpd = TRUE)
# #   segments(block_start, 0, block_start, 0.1, col = "blue", lwd = 2)
# #   segments(block_end,   0, block_end,   0.1, col = "blue", lwd = 2)
# #   abline(h = 0.1, col = "blue", lty = 2, lwd = 1.5)
# #   duration_days <- as.integer(block_end - block_start) + 1
# #   text(mean(c(block_start, block_end)), 0.05,
# #        paste0(duration_days, " Days"), col = "blue", cex = 0.9, font = 2)
# # }
# # 
# # idx <- which(classified_data_us$share >= 0.01)[1]
# # if (!is.na(idx)) {
# #   x_pt <- classified_data_us$date[idx]
# #   y_pt <- classified_data_us$share[idx]
# # 
# #   points(x_pt, y_pt, col = "red", pch = 16, cex = 1.2)
# #   arrows(x_pt, y_pt + 0.12, x_pt, y_pt + 0.01, length = 0.1, col = "red", lwd = 1.5)
# #   text(x_pt, y_pt + 0.14, expression("1"^{st}~"over 1% share"),
# #        col = "red", cex = 0.9, pos = 3)
# # }






# # Sep 26, 2025
# lin_name <- unique(classified_data_us$lineage)
# 
# idx <- which(classified_data_us$share >= 0.01)[1]
# 
# subset_list <- list()
# if (!is.na(idx)) {
#   first_date <- classified_data_us$date[idx]
#   
#   weeks <- c(2, 3, 4, 5)
#   
#   subset_list <- lapply(weeks, function(w) {
#     classified_data_us %>%
#       filter(date >= first_date & date < first_date + 7 * w) %>%
#       mutate(window = paste0(w, "_weeks"))
#   })
#   
#   names(subset_list) <- paste0(weeks, "_weeks")
# }
# 
# lin_name <- unique(classified_data_us$lineage)
# cty_name <- unique(classified_data_us$country)
# 
# idx <- which(classified_data_us$share >= 0.01)[1]
# 
# if (!is.na(idx)) {
#   first_date <- classified_data_us$date[idx]
#   weeks <- c(2, 3, 4, 5)
#   
#   subset_list <- lapply(weeks, function(w) {
#     classified_data_us %>%
#       filter(date >= first_date & date < first_date + 7 * w) %>%
#       mutate(window = paste0(w, "_weeks"))
#   })
#   names(subset_list) <- paste0(weeks, "_weeks")
#   
#   oldpar <- par(no.readonly = TRUE)
#   on.exit(par(oldpar))
#   par(mfrow = c(4, 1), mar = c(4, 4, 3, 1))
#   
#   for (i in seq_along(subset_list)) {
#     df_sub <- subset_list[[i]]
#     w <- weeks[i]
#     
#     plot(df_sub$date, df_sub$share,
#          type = "o", ylim = c(0, 0.1),
#          xlab = "", ylab = "Share (%)",
#          pch = 16, cex = 0.6,
#          xaxt = "n",
#          main = paste0(w, " weeks"),
#          yaxt = "n",
#          xlim = c(first_date, first_date + 7*4))
#     
#     axis(2, at = seq(0, 0.1, by = 0.02),
#          labels = paste0(seq(0, 10, by = 2), "%"))
#     
#     axis.Date(1, at = seq(first_date, first_date + 7*4, by = "1 week"),
#               format = "%b %d")
#     
#     points(df_sub$date[1], df_sub$share[1], col = "red", pch = 16, cex = 1.5)
#   }
# }


# ------------------------------------------------------------------------------
# Build early-period aligned wide tables for multiple input windows
# ------------------------------------------------------------------------------
input_days_vec <- c(14, 21, 28, 35)
res_list <- setNames(
  lapply(input_days_vec, function(d) {
    prepare_model_data(
      data               = classified_data,
      input_days         = d,
      dominant_threshold = 0.5,
      min_share          = 0.01,
      values_fill        = NA_real_
    )
  }),
  paste0("input_", input_days_vec)
)

feat_list <- lapply(c(14, 21, 28, 35), function(d) {
  extract_features_from_model_data(res_list[[paste0("input_", d)]]$X, add_auto_features = TRUE)
})
names(feat_list) <- paste0("feat_", c(14, 21, 28, 35))
# saveRDS(feat_list, file = "code2/feat_list.rds")

# feat_list$feat_14 %>% dim()
# feat_list$feat_21 %>% dim()
# feat_list$feat_28 %>% dim()
# feat_list$feat_35 %>% dim()

# # feat_list$feat_21 %>% colnames()
# [1] "country"                                         "lineage"                                        
# [3] "slope_log"                                       "slope_y_1_7"                                    
# [5] "peak_val"                                        "auc_norm"                                       
# [7] "longest_inc_run"                                 "zero_count"                                     
# [9] "c22_DN_HistogramMode_5"                          "c22_DN_HistogramMode_10"                        
# [11] "c22_CO_f1ecac"                                   "c22_CO_FirstMin_ac"                             
# [13] "c22_CO_HistogramAMI_even_2_5"                    "c22_CO_trev_1_num"                              
# [15] "c22_MD_hrv_classic_pnn40"                        "c22_SB_BinaryStats_mean_longstretch1"           
# [17] "c22_SB_TransitionMatrix_3ac_sumdiagcov"          "c22_PD_PeriodicityWang_th0_01"                  
# [19] "c22_CO_Embed2_Dist_tau_d_expfit_meandiff"        "c22_IN_AutoMutualInfoStats_40_gaussian_fmmi"    
# [21] "c22_FC_LocalSimple_mean1_tauresrat"              "c22_DN_OutlierInclude_p_001_mdrmd"              
# [23] "c22_DN_OutlierInclude_n_001_mdrmd"               "c22_SP_Summaries_welch_rect_area_5_1"           
# [25] "c22_SB_BinaryStats_diff_longstretch0"            "c22_SB_MotifThree_quantile_hh"                  
# [27] "c22_SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1" "c22_SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1"     
# [29] "c22_SP_Summaries_welch_rect_centroid"            "c22_FC_LocalSimple_mean3_stderr" 


# # Silhouette + Jenks to select the number of groups 
# vol_groups <- group_countries_by_volatility_index(
#   filtered_data = filtered_data,
#   n_groups      = 4,          
#   label_sep     = "-",        
#   breaks_method = "jenks",    
#   seed          = 123
# )
# 
# vol_groups_view <- vol_groups %>%
#   arrange(desc(volatility_index)) %>%
#   dplyr::select(country, volatility_index, volatility_group)
# # print(vol_groups_view, n = Inf)
# 
# saveRDS(vol_groups_view, file = "code/vol_groups_view.rds")


# df_1 <- res_list$input_14$combo_flags %>% 
#   filter(window_len_ge_input == FALSE | 
#          reached_dominant_in_window == TRUE | 
#          any_above_min_share == FALSE) %>% 
#   select(country, lineage, 
#          window_len_ge_input, reached_dominant_in_window, any_above_min_share) %>% 
#   mutate(input_14 = "yes")
# 
# df_2 <- res_list$input_21$combo_flags %>% 
#   filter(window_len_ge_input == FALSE | 
#            reached_dominant_in_window == TRUE | 
#            any_above_min_share == FALSE) %>% 
#   select(country, lineage, 
#          window_len_ge_input, reached_dominant_in_window, any_above_min_share) %>% 
#   mutate(input_21 = "yes")
# 
# df_3 <- res_list$input_28$combo_flags %>% 
#   filter(window_len_ge_input == FALSE | 
#            reached_dominant_in_window == TRUE | 
#            any_above_min_share == FALSE) %>% 
#   select(country, lineage, 
#          window_len_ge_input, reached_dominant_in_window, any_above_min_share) %>% 
#   mutate(input_28 = "yes")
# 
# final_df <- res_list$input_35$combo_flags %>% 
#   filter(window_len_ge_input == FALSE | 
#            reached_dominant_in_window == TRUE | 
#            any_above_min_share == FALSE) %>% 
#   select(country, lineage, first_share_date,
#          window_len_ge_input, reached_dominant_in_window, any_above_min_share) %>% 
#   mutate(exclusion = case_when(
#     paste(country, lineage) %in% paste(df_1$country, df_1$lineage) ~ "input_14",
#     paste(country, lineage) %in% paste(df_2$country, df_2$lineage) ~ "input_21",
#     paste(country, lineage) %in% paste(df_3$country, df_3$lineage) ~ "input_28",
#     TRUE ~ "input_35"
#   ))
# 
# final_df %>% 
#   filter(!(paste(country, lineage) %in% paste(df_3$country, df_3$lineage)))



