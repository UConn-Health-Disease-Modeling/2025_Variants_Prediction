# clear the environment
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(caret)
  library(Cubist)
  library(iml)
  library(ggplot2)
})

set.seed(2025)

# ---------- Load & align data (use 35-day features) ----------
feat_list    <- readRDS("code/feat_list.rds")
measurements <- readRDS("code/all_data.rds")[["measurements"]] %>%
  dplyr::select(country, lineage, peak_share)

split_one <- function(feat_df, meas_df) {
  df <- feat_df %>% inner_join(meas_df, by = c("country","lineage"))
  id_cols <- c("country","lineage"); y_col <- "peak_share"
  feat_cols <- setdiff(names(df), c(id_cols, y_col))
  df[feat_cols] <- lapply(df[feat_cols], as.numeric)
  all_na <- apply(df[feat_cols], 1, function(r) all(is.na(r)))
  if (any(all_na)) df <- df[!all_na, , drop = FALSE]
  list(
    X = as.matrix(df[, feat_cols, drop = FALSE]),
    Y = as.numeric(df[[y_col]]),
    ids = df[, id_cols, drop = FALSE],
    feature_names = feat_cols
  )
}

sp <- split_one(feat_list$feat_28, measurements)

# ---------- Median imputation ----------
impute_median <- function(X, med = NULL) {
  if (is.null(med)) {
    med <- apply(X, 2, function(x) {
      m <- suppressWarnings(stats::median(x, na.rm = TRUE))
      if (!is.finite(m)) NA_real_ else m
    })
  }
  Xi <- X
  for (j in seq_along(med)) {
    idx <- is.na(Xi[, j])
    if (any(idx)) Xi[idx, j] <- med[j]
  }
  list(X = Xi, med = med)
}

imp   <- impute_median(sp$X)
X_imp <- imp$X
Y     <- sp$Y

# ---------- Tune & fit Cubist ----------
ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE, allowParallel = TRUE)
grid <- expand.grid(committees = c(10, 50, 100), neighbors = c(0, 3, 5))

tune_fit <- caret::train(
  x = as.data.frame(X_imp), y = Y,
  method = "cubist",
  trControl = ctrl,
  tuneGrid  = grid,
  metric    = "RMSE"
)

bp <- tune_fit$bestTune
final_model <- Cubist::cubist(
  x = as.data.frame(X_imp),
  y = Y,
  committees = bp$committees,
  neighbors  = bp$neighbors
)

# ---------- SHAP via iml (Shapley) ----------
X_df <- as.data.frame(X_imp)
pred_fun_iml <- function(m, newdata) {  # prediction wrapper for iml
  nd <- as.data.frame(newdata)
  nd[] <- lapply(nd, as.numeric)
  as.numeric(predict(m, newdata = nd))
}

predictor <- Predictor$new(
  model = final_model,
  data  = X_df,
  y     = Y,
  predict.function = pred_fun_iml,
  type  = "response"
)

set.seed(2025)
n_explain   <- min(200, nrow(X_df))           # subsample for speed
idx_explain <- sample.int(nrow(X_df), n_explain)

get_one_shap <- function(i) {
  sh <- Shapley$new(predictor, x.interest = X_df[i, , drop = FALSE], sample.size = 100)
  out <- sh$results
  out$case_id <- i
  out
}

shap_all <- bind_rows(lapply(idx_explain, get_one_shap))
top_feats <- shap_all %>%
  group_by(feature) %>%
  summarise(mean_abs_shap = mean(abs(phi), na.rm = TRUE)) %>%  # 改成 mean_abs_shap
  arrange(desc(mean_abs_shap)) %>%
  slice_head(n = 10)

feature_labels <- c(
  peak_val = "Early peak value",
  zero_count = "non-zero-counts",
  c22_SB_MotifThree_quantile_hh = "Motif frequency (3-state)",
  c22_DN_OutlierInclude_p_001_mdrmd = "Outlier prevalence",
  c22_FC_LocalSimple_mean3_stderr = "Local fluctuation (stderr)",
  c22_SP_Summaries_welch_rect_area_5_1 = "Spectral Welch area",
  c22_SB_BinaryStats_mean_longstretch1 = "Longest run (binary)",
  c22_CO_HistogramAMI_even_2_5 = "Auto mutual information (hist)",
  c22_PD_PeriodicityWang_th0_01 = "Periodicity score",
  c22_SB_TransitionMatrix_3ac_sumdiagcov = "Transition-matrix covariance"
)

# feature_labels <- c(
#   peak_val = "Early peak value",
#   zero_count = "Non-zero counts",
#   auc_norm = "Normalized AUC",
#   c22_IN_AutoMutualInfoStats_40_gaussian_fmmi = "Auto mutual info (lag=40, Gaussian)",
#   c22_DN_OutlierInclude_p_001_mdrmd = "Outlier prevalence",
#   c22_SB_BinaryStats_mean_longstretch1 = "Longest run (binary)",
#   c22_CO_trev_1_num = "Time reversal asymmetry (lag=1)",
#   c22_SB_TransitionMatrix_3ac_sumdiagcov = "Transition-matrix covariance",
#   c22_SB_MotifThree_quantile_hh = "Motif frequency (3-state)",
#   c22_CO_Embed2_Dist_tau_d_expfit_meandiff = "Embedding distance (exp fit, mean diff)"
# )

top_feats <- top_feats %>%
  mutate(feature_clean = recode(feature, !!!feature_labels))

# ---------- Global importance (mean |SHAP|) & top-10 ----------
p_imp <- ggplot(top_feats, aes(x = reorder(feature_clean, mean_abs_shap), 
                               y = mean_abs_shap)) +
  geom_col(width = 0.75, fill = "#E31A1C") +
  coord_flip() +
  expand_limits(y = max(top_feats$mean_abs_shap) * 1.15) +
  labs(x = NULL, y = "Mean |SHAP|", title = "Top 10 Global Importance (mean |SHAP|)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_imp)

# ---------- Beeswarm (top-10, within-feature color scaling) ----------
if (!"x.recoded" %in% names(shap_all)) {
  shap_all <- shap_all %>%
    left_join(
      X_df %>%
        mutate(case_id = row_number()) %>%
        pivot_longer(-case_id, names_to = "feature", values_to = "x.recoded"),
      by = c("case_id","feature")
    )
}

shap_all <- shap_all %>%
  filter(feature %in% top_feats$feature) %>%
  mutate(
    feature_clean = recode(feature, !!!feature_labels),
    feature_clean = factor(feature_clean, levels = rev(top_feats$feature_clean))
  ) %>%
  group_by(feature_clean) %>%
  mutate(
    val_pct  = percent_rank(x.recoded),
    col_val  = val_pct^0.6,
    zero_var = sd(x.recoded, na.rm = TRUE) == 0
  ) %>%
  ungroup()

shap_all$col_val[shap_all$zero_var] <- NA

p_bee <- ggplot(shap_all, aes(x = feature_clean, y = phi,
                              color = col_val, alpha = col_val, size = col_val)) +
  geom_jitter(width = 0.22, height = 0, shape = 16, stroke = 0) +
  scale_color_gradientn(
    colors = c("#2166AC", "#6EA6CD", "#F7A789", "#B2182B"),
    values = c(0, 0.5, 0.8, 1),
    na.value = "grey70",
    breaks = c(0, 1),
    labels = c("Low", "High"),
    name = "Feature\nvalue"
  ) +
  scale_alpha(range = c(0.25, 0.95), guide = "none") +
  scale_size(range = c(0.6, 1.8), guide = "none") +
  labs(x = NULL, y = "SHAP value", title = "SHAP Beeswarm (Top 10 Features)") +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.key.width = unit(0.4, "cm"),
    legend.title = element_text(angle = 0, vjust = 0.5, hjust = 0.5)
  )

print(p_bee)

