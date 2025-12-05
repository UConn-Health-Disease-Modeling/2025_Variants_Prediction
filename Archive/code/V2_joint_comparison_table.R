rm(list = ls())
library(dplyr)


res_peak     <- readRDS("code/peak_results_1006.rds")$feat_28_ensemble$result_df
res_duration <- readRDS("code/duration_results_1006.rds")$feat_21_ensemble$final_df
all_data     <- readRDS("code/all_data.rds")
measurements      <- all_data$measurements

joint_compare <- measurements %>%
  dplyr::select(country, lineage, start_date, end_date) %>%
  left_join(res_peak %>% dplyr::select(country, lineage, peak_share_pred, peak_share), 
            by = c("country", "lineage")) %>%
  left_join(res_duration %>% dplyr::select(country, lineage, predicted_days_above_10_cat, true_days_above_10_cat),
            by = c("country", "lineage")) %>%
  tidyr::drop_na()

mean(abs(joint_compare$peak_share - joint_compare$peak_share_pred))

joint_compare_filtered <- joint_compare %>%
  dplyr::filter(
    peak_share > 0.3 & true_days_above_10_cat == "100+"
  ) %>%
  dplyr::arrange(dplyr::desc(country), start_date)


joint_compare %>%
  dplyr::filter(
    (peak_share > 0.3 & true_days_above_10_cat == "100+")
  ) %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(
    n_variants = dplyr::n(),
    earliest_start = min(start_date, na.rm = TRUE),
    latest_end = max(end_date, na.rm = TRUE),
    .groups = "drop"
  )

#    country        n_variants earliest_start latest_end
#    <chr>               <int> <date>         <date>    
#  1 Australia              32 2021-06-16     2024-05-03
#  2 Austria                18 2021-06-14     2023-12-01
#  3 Brazil                 27 2020-06-01     2024-02-21
#  4 Canada                 76 2020-03-16     2024-09-01
#  5 Denmark                54 2021-06-21     2023-04-23
#  6 France                 61 2020-03-02     2024-07-02
#  7 Germany                86 2020-03-18     2024-01-17
#  8 India                  21 2020-04-27     2024-01-08
#  9 Italy                  27 2020-07-20     2024-01-15
# 10 Japan                  75 2020-03-23     2024-08-17
# 11 South Korea            25 2021-02-09     2024-03-27
# 12 Spain                  38 2020-10-29     2024-08-11
# 13 Sweden                 36 2020-05-12     2024-07-19
# 14 United Kingdom         63 2020-03-04     2024-08-25
# 15 United States          91 2020-03-02     2024-09-01

library(openxlsx)

feature_df <- data.frame(
  "#" = 1:24,
  "Feature name" = c(
    "DN_HistogramMode_5","DN_HistogramMode_10","DN_OutlierInclude_p_001_mdrmd","DN_OutlierInclude_n_001_mdrmd",
    "first1e_acf_tau","firstMin_acf","SP_Summaries_welch_rect_area_5_1","SP_Summaries_welch_rect_centroid",
    "FC_LocalSimple_mean3_stderr","FC_LocalSimple_mean1_tauresrat","MD_hrv_classic_pnn40","SB_BinaryStats_mean_longstretch1",
    "SB_BinaryStats_diff_longstretch0","SB_MotifThree_quantile_hh","CO_HistogramAMI_even_2_5","CO_trev_1_num",
    "IN_AutoMutualInfoStats_40_gaussian_fmmi","SB_TransitionMatrix_3ac_sumdiagcov","PD_PeriodicityWang_th001",
    "CO_Embed2_Dist_tau_d_expfit_meandiff","SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1","SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1",
    "DN_Mean","DN_Spread_Std"
  ),
  "Short name" = c(
    "mode_5","mode_10","outlier_timing_pos","outlier_timing_neg",
    "acf_timescale","acf_first_min","low_freq_power","centroid_freq",
    "forecast_error","whiten_timescale","high_fluctuation","stretch_high",
    "stretch_decreasing","entropy_pairs","ami2","trev",
    "ami_timescale","transition_variance","periodicity",
    "embedding_dist","rs_range","dfa","mean","std"
  ),
  "Category" = c(
    "Distribution shape","Distribution shape","Extreme event timing","Extreme event timing",
    "Linear autocorrelation","Linear autocorrelation","Linear autocorrelation","Linear autocorrelation",
    "Simple forecasting","Incremental differences","Incremental differences","Symbolic",
    "Symbolic","Symbolic","Nonlinear autocorrelation","Nonlinear autocorrelation",
    "Linear autocorrelation","Symbolic","Linear autocorrelation",
    "Linear autocorrelation","Self-affine scaling","Self-affine scaling","Catch24 extension","Catch24 extension"
  ),
  "Description" = c(
    "5-bin histogram mode","10-bin histogram mode","Positive outlier timing","Negative outlier timing",
    "First 1/e crossing of the ACF","First minimum of the ACF","Power in lowest 20% frequencies","Centroid frequency",
    "Error of 3-point rolling mean forecast","Change in autocorr. timescale after differencing","Proportion of high incremental changes","Longest stretch of above-mean values",
    "Longest stretch of decreasing values","Entropy of successive pairs in symbolized series","Histogram-based AMI (lag 2, 5 bins)","Time reversibility",
    "First minimum of AMI function","Transition matrix column variance","Wang's periodicity metric",
    "Exponential fit to embedding dist. distribution","Rescaled range fluctuation analysis (low-scale)","Detrended fluctuation analysis (low-scale)","Mean of raw series","Standard deviation of raw series"
  ),
  stringsAsFactors = FALSE
)
