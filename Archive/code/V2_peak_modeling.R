# clear the environment
rm(list = ls())

name_map <- c(
  slope_log  = "log_slope",
  slope_y_1_7 = "short_slopes",
  longest_inc_run = "inc_run",
  zero_count = "nonzero_days",
  
  c22_DN_HistogramMode_5 = "mode_5",
  c22_DN_HistogramMode_10 = "mode_10",
  c22_CO_f1ecac = "acf_timescale",
  c22_CO_FirstMin_ac = "acf_first_min",
  c22_CO_HistogramAMI_even_2_5 = "ami2",
  c22_CO_trev_1_num = "trev",
  c22_MD_hrv_classic_pnn40 = "high_fluctuation",
  c22_SB_BinaryStats_mean_longstretch1 = "stretch_high",
  c22_SB_TransitionMatrix_3ac_sumdiagcov = "transition_variance",
  c22_PD_PeriodicityWang_th0_01 = "periodicity",
  c22_CO_Embed2_Dist_tau_d_expfit_meandiff = "embedding_dist",
  c22_IN_AutoMutualInfoStats_40_gaussian_fmmi = "ami_timescale",
  c22_FC_LocalSimple_mean1_tauresrat = "whiten_timescale",
  c22_DN_OutlierInclude_p_001_mdrmd = "outlier_timing_pos",
  c22_DN_OutlierInclude_n_001_mdrmd = "outlier_timing_neg",
  c22_SP_Summaries_welch_rect_area_5_1 = "low_freq_power",
  c22_SB_BinaryStats_diff_longstretch0 = "stretch_decreasing",
  c22_SB_MotifThree_quantile_hh = "entropy_pairs",
  c22_SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1 = "rs_range",
  c22_SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1 = "dfa",
  c22_SP_Summaries_welch_rect_centroid = "centroid_freq",
  c22_FC_LocalSimple_mean3_stderr = "forecast_error"
)

inv_logit_peak <- function(peak_share_logit, n) {
  eps <- 0.5 / (n + 1)
  peak_share_est <- plogis(peak_share_logit)  
  peak_share <- pmin(pmax(peak_share_est, eps), 1 - eps)
  return(peak_share)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(Cubist)
  library(SuperLearner)
  library(glmnet)        # for SL.glmnet
  library(randomForest)  # for SL.randomForest
  library(xgboost)       # for SL.xgboost
  library(MASS)          # for SL.bayesglm
  library(gam)           # for SL.gam
  library(kknn)          # for SL.knn
  library(e1071)         # for SL.svm
  library(arm)
  library(tidyr)
  library(rpart)
})

feat_list    <- readRDS("code/feat_list.rds")
measurements <- readRDS("code/all_data.rds")[["measurements"]] %>%
  dplyr::select(country, lineage, peak_share, peak_share_logit)

run_superlearner <- function(feat_data,
                             measurements,
                             learners = c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.xgboost"),
                             seed = 123) {
  library(dplyr)
  library(SuperLearner)
  
  set.seed(seed)
  
  model_df <- feat_data %>%
    dplyr::select(-auc_norm, -peak_val) %>%
    left_join(measurements, by = c("country", "lineage")) %>%
    na.omit()
  
  train_df <- model_df %>%
    dplyr::select(-country, -lineage)
  
  Y <- train_df$peak_share_logit
  X <- train_df %>% dplyr::select(-peak_share_logit, -peak_share)
  
  sl_fit <- SuperLearner(
    Y = Y,
    X = X,
    SL.library = learners,
    family = gaussian(),
    cvControl = list(V = 5)
  )
  
  n <- nrow(X)
  inv_logit_peak <- function(peak_share_logit, n) {
    eps <- 0.5 / (n + 1)
    peak_share_est <- plogis(peak_share_logit)
    peak_share <- pmin(pmax(peak_share_est, eps), 1 - eps)
    return(peak_share)
  }
  
  model_df$peak_share_pred <- inv_logit_peak(sl_fit$SL.predict, n)
  result_df <- model_df %>%
    dplyr::select(country, lineage, peak_share, peak_share_pred)
  
  rmse <- sqrt(mean((result_df$peak_share - result_df$peak_share_pred)^2))
  ss_res <- sum((result_df$peak_share - result_df$peak_share_pred)^2)
  ss_tot <- sum((result_df$peak_share - mean(result_df$peak_share))^2)
  r2 <- 1 - ss_res / ss_tot
  
  list(
    sl_fit = sl_fit,
    result_df = result_df,
    rmse = rmse,
    r2 = r2,
    coef = sl_fit$coef
  )
}


# ------------------- Setup -------------------
feature_sets <- c("feat_14", "feat_21", "feat_28", "feat_35")

learner_combos <- list(
  "glm" = "SL.glm",
  "bayesglm" = "SL.bayesglm",
  "gam" = "SL.gam",
  "elasticnet" = "SL.glmnet",
  "svm" = "SL.svm",
  "nnet" = "SL.nnet",
  "CART" = "SL.rpart",
  "ensemble" = c("SL.randomForest", "SL.glmnet", "SL.xgboost")
)

results_list <- list()

# ------------------- Loop -------------------
for (fname in feature_sets) {
  for (mname in names(learner_combos)) {
    cat("Running SuperLearner for:", fname, "-", mname, "\n")
    
    res <- run_superlearner(
      feat_data = feat_list[[fname]],
      measurements = measurements,
      learners = learner_combos[[mname]],
      seed = 2025
    )
    
    results_list[[paste(fname, mname, sep = "_")]] <- res
  }
}

# saveRDS(results_list, "peak_results_1006.rds")

# ------------------- Summary -------------------
summary_df <- data.frame(
  feature_set = sub("^(feat_\\d+)_.*", "\\1", names(results_list)),
  model = sub("^feat_\\d+_(.*)$", "\\1", names(results_list)),
  RMSE = sapply(results_list, function(x) x$rmse),
  R2 = sapply(results_list, function(x) x$r2)
)


# summary_wide <- summary_df %>% mutate(model = factor(model, levels = c("glm", "bayesglm", "gam", "elasticnet", "svm", "nnet", "CART", "ensemble"))) %>%
#   arrange(model, feature_set) %>%
#   pivot_wider(
#     names_from = feature_set,
#     values_from = c(RMSE, R2)
#   ) %>%
#   dplyr::select(model,
#                 RMSE_feat_14, R2_feat_14,
#                 RMSE_feat_21, R2_feat_21,
#                 RMSE_feat_28, R2_feat_28,
#                 RMSE_feat_35, R2_feat_35)


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(scales)
})

df_list <- lapply(names(results_list), function(nm) {
  if (grepl("^feat_28", nm)) {
    tmp <- results_list[[nm]]$result_df
    tmp$model <- gsub("feat_28_", "", nm)
    tmp
  } else NULL
})
df_list <- df_list[!sapply(df_list, is.null)]
all_df <- bind_rows(df_list)

country_rmse <- all_df %>%
  group_by(model, country) %>%
  summarise(rmse = sqrt(mean((peak_share - peak_share_pred[,1])^2, na.rm = TRUE)), .groups = "drop")

model_name_map <- c(
  "glm" = "GLM",
  "bayesglm" = "BayesGLM",
  "gam" = "GAM",
  "elasticnet" = "Elastic Net",
  "svm" = "SVM",
  "nnet" = "Neural Net",
  "CART" = "CART",
  "ensemble" = "SuperLearner"
)

order_models <- c("GLM", "BayesGLM", "GAM", "Elastic Net",
                  "SVM", "Neural Net", "CART", "SuperLearner")

heat_data <- country_rmse %>%
  mutate(model = model_name_map[model]) %>%
  mutate(model = factor(model, levels = order_models))

ggplot(heat_data, aes(x = model, y = reorder(country, desc(country)), fill = rmse)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rmse)), size = 3, color = "black") +
  scale_fill_gradient(low = "#8b0000", high = "#fbe6d4", limits = c(0.1, 0.27), oob = squish) +
  labs(title = "RMSE by country and model (28 Days Input)", fill = "RMSE") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )









# library(SuperLearner)
# library(iml)
# library(ggplot2)
# library(DALEX)
# pred_fun <- function(model, newdata) {
#   predict(model, newdata = newdata)$pred
# }
# 
# explainer_sl <- explain(
#   model = sl_fit,
#   data = X,
#   y = Y,
#   predict_function = pred_fun,
#   label = "SuperLearner"
# )
# 
# mp <- model_parts(explainer_sl)
# plot(mp, show_boxplots = FALSE) + ggtitle("Feature Importance", "")
# 
# 
# bd <- predict_parts(explainer_sl, new_observation = X[213, ])
# 
# labs <- bd$variable
# has_eq <- grepl("=", labs)
# vals <- regmatches(labs[has_eq], regexpr("=\\s*-?\\d+\\.?\\d*", labs[has_eq]))
# nums <- as.numeric(sub("=\\s*", "", vals))
# nums_fmt <- formatC(nums, format = "f", digits = 2, drop0trailing = TRUE)
# labs[has_eq] <- mapply(function(orig, newval) {
#   sub("=\\s*-?\\d+\\.?\\d*", paste0("= ", newval), orig)
# }, labs[has_eq], nums_fmt, USE.NAMES = FALSE)
# 
# bd$variable <- labs
# 
# plot(bd)
















# example model
library(xgboost)
library(dplyr)

feat_df <- feat_list$feat_35 %>%
  dplyr::select(-auc_norm, -peak_val) %>%
  left_join(measurements, by = c("country", "lineage")) %>%
  na.omit()

Y <- feat_df$peak_share_logit
X <- feat_df %>%
  dplyr::select(-country, -lineage, -peak_share_logit, -peak_share) %>%
  as.matrix()

dtrain <- xgb.DMatrix(data = X, label = Y)

set.seed(123)
xgb_model <- xgboost(
  data = dtrain,
  objective = "reg:squarederror",
  nrounds = 300,
  max_depth = 4,
  eta = 0.05,
  subsample = 0.8,
  colsample_bytree = 0.8,
  verbose = 0
)

importance_matrix <- xgb.importance(model = xgb_model)
importance_df <- xgb.importance(model = xgb_model) %>%
  as_tibble() %>%
  dplyr::mutate(
    Feature = ifelse(
      Feature %in% names(name_map),
      name_map[Feature],
      Feature
    )
  ) %>%
  dplyr::arrange(desc(Gain)) %>%
  dplyr::slice_head(n = 10)

library(ggplot2)
ggplot(importance_df, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_col(fill = "#E31A1C", width = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Gain", title = "Top 10 Feature Importance") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())







library(iml)
library(ggplot2)
library(dplyr)

pred_fun <- function(model, newdata) {
  as.numeric(predict(model, newdata = as.matrix(newdata)))
}

predictor <- Predictor$new(
  model = xgb_model,
  data = as.data.frame(X),
  y = Y,
  predict.function = pred_fun
)

set.seed(123)
n_explain <- min(200, nrow(X))
idx_explain <- sample.int(nrow(X), n_explain)

get_one_shap <- function(i) {
  sh <- Shapley$new(predictor, x.interest = as.data.frame(X[i, , drop = FALSE]), sample.size = 100)
  out <- sh$results
  out$case_id <- i
  out
}

shap_all <- bind_rows(lapply(idx_explain, get_one_shap))

shap_all <- shap_all %>%
  left_join(
    as.data.frame(X) %>%
      mutate(case_id = row_number()) %>%
      tidyr::pivot_longer(-case_id, names_to = "feature", values_to = "x.recoded"),
    by = c("case_id","feature")
  ) %>%
  mutate(
    feature = ifelse(feature %in% names(name_map), name_map[feature], feature)
  )

# 使用 importance_df 中的特征
selected_feats <- importance_df$Feature

shap_all <- shap_all %>%
  filter(feature %in% selected_feats) %>%
  mutate(feature = factor(feature, levels = rev(selected_feats))) %>%
  group_by(feature) %>%
  mutate(
    val_scaled = (x.recoded - median(x.recoded, na.rm = TRUE)) / sd(x.recoded, na.rm = TRUE),
    val_pct = pnorm(val_scaled),
    col_val = val_pct
  ) %>%
  ungroup() 

ggplot(shap_all, aes(x = feature, y = phi, color = col_val)) +
  geom_jitter(width = 0.25, height = 0, size = 1, alpha = 0.8) +
  scale_color_gradientn(
    colors = c("#2166AC", "#6EA6CD", "#F7A789", "#B2182B"),
    values = c(0, 0.3, 0.5, 1),
    name = "Feature value",
    breaks = c(0, 1),
    labels = c("Lower", "Higher"),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(5, "cm"),
      barwidth = unit(0.5, "cm")
    )
  ) +
  coord_flip() +
  labs(x = NULL, y = "SHAP value", title = "SHAP Beeswarm (Top 10 XGBoost Features)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )