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
  select(country, lineage, days_above_10_cat)

split_one <- function(feat_df, meas_df,
                      y_col = "days_above_10_cat",
                      id_cols = c("country","lineage"),
                      ordered_levels = c("0","1–30","31–100","100+"),
                      make_ordinal = FALSE) {
  df <- dplyr::inner_join(feat_df, meas_df, by = id_cols)
  feat_cols <- setdiff(names(df), c(id_cols, y_col))
  df[feat_cols] <- lapply(df[feat_cols], function(x) if (is.numeric(x)) x else suppressWarnings(as.numeric(x)))
  all_na_row <- apply(df[feat_cols], 1, function(r) all(is.na(r)))
  if (any(all_na_row)) df <- df[!all_na_row, , drop = FALSE]
  y_raw <- df[[y_col]]
  if (is.numeric(y_raw)) y_raw <- as.character(y_raw)
  observed <- unique(as.character(y_raw))
  use_levels <- if (length(intersect(ordered_levels, observed)) >= 2) ordered_levels[ordered_levels %in% observed] else sort(observed)
  y_factor <- factor(as.character(y_raw), levels = use_levels, ordered = make_ordinal)
  y_index_1based <- as.integer(y_factor)
  y_index_0based <- y_index_1based - 1L
  X <- as.matrix(df[, feat_cols, drop = FALSE])
  list(
    X = X,
    Y_factor = y_factor,
    Y_index_1based = y_index_1based,
    Y_index_0based = y_index_0based,
    y_levels = levels(y_factor),
    ids = df[, id_cols, drop = FALSE],
    feature_names = feat_cols,
    is_ordered = is.ordered(y_factor)
  )
}

days <- c(28, 35)
available <- intersect(paste0("feat_", days), names(feat_list))
splits <- setNames(
  lapply(available, function(nm) split_one(feat_list[[nm]], measurements)),
  sub("^feat_", "split_", available)
)

# splits1_countries <- c(
#   "United Kingdom","United States",
#   "Germany","Japan","Canada","France"
# )
splits2_countries <- c(
  "Denmark","Sweden","India","Brazil",
  "Australia","Spain","Italy","Austria","South Korea", 
  "United Kingdom","United States",
  "Germany","Japan","Canada","France"
)

make_country_splits <- function(splits, countries) {
  lapply(splits, function(sp) {
    keep <- sp$ids$country %in% countries
    lapply(sp, function(x) {
      if (is.data.frame(x) || is.matrix(x)) x[keep, , drop = FALSE] else x[keep]
    })
  })
}

# splits1 <- make_country_splits(splits, splits1_countries)
splits2 <- make_country_splits(splits, splits2_countries)

macro_metrics <- function(truth, pred) {
  truth <- factor(truth)
  pred  <- factor(pred, levels = levels(truth))
  cm <- table(truth, pred)
  tp <- diag(cm)
  fp <- colSums(cm) - tp
  fn <- rowSums(cm) - tp
  precision <- tp / pmax(tp + fp, 1)
  recall <- tp / pmax(tp + fn, 1)
  f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
  acc <- sum(tp) / sum(cm)
  bal_acc <- mean(recall, na.rm = TRUE)
  c(accuracy = acc, macro_f1 = mean(f1, na.rm = TRUE), balanced_accuracy = bal_acc)
}

fit_and_oof_elasticnet <- function(sp, seed = 1L) {
  X <- as.data.frame(sp$X)
  Y <- if (!is.null(sp$Y_factor)) sp$Y_factor else factor(sp$Y)
  if (!is.factor(Y)) Y <- factor(Y)
  if (nlevels(Y) < 2) stop("Y has fewer than 2 classes.")
  
  keep_cols <- colSums(!is.na(X)) > 0
  X <- X[, keep_cols, drop = FALSE]
  
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = FALSE,
    summaryFunction = defaultSummary,
    savePredictions = "final",
    allowParallel = TRUE
  )
  
  grid <- expand.grid(
    alpha  = seq(0, 1, by = 0.25),
    lambda = 10^seq(0, -4, length.out = 30)
  )
  
  set.seed(seed)
  en_fit <- caret::train(
    x = X,
    y = Y,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid  = grid,
    metric    = "Accuracy",
    preProcess = c("medianImpute","zv","center","scale"),
    family = "multinomial"
  )
  
  best_params <- en_fit$bestTune
  
  oof <- en_fit$pred
  oof <- oof[oof$alpha == best_params$alpha & 
               abs(oof$lambda - best_params$lambda) < .Machine$double.eps^0.5, , drop = FALSE]
  oof <- oof[, c("rowIndex","Resample","obs","pred")]
  colnames(oof) <- c("rowIndex","fold","y_true","y_pred")
  oof <- oof[order(oof$rowIndex), ]
  
  ids_aligned <- sp$ids[oof$rowIndex, , drop = FALSE]
  pred_df <- dplyr::bind_cols(ids_aligned, oof[, c("fold","y_true","y_pred")])
  
  mets <- macro_metrics(truth = pred_df$y_true, pred = pred_df$y_pred)
  
  list(
    best_params       = best_params,
    accuracy          = unname(mets["accuracy"]),
    macro_f1          = unname(mets["macro_f1"]),
    balanced_accuracy = unname(mets["balanced_accuracy"]),
    oof_predictions   = pred_df,
    varimp            = caret::varImp(en_fit),
    final_model       = en_fit$finalModel,  
    caret_fit         = en_fit
  )
}

results1 <- lapply(splits1, fit_and_oof_elasticnet, seed = 1)




sp        <- splits1$split_28
res_split <- results1$split_28
final_model <- res_split$caret_fit

X_df <- as.data.frame(sp$X)
X_df[] <- lapply(X_df, as.numeric)  
Y    <- sp$Y_factor

pred_fun_iml <- function(m, newdata) {
  nd <- as.data.frame(newdata)
  nd[] <- lapply(nd, as.numeric)
  predict(m, newdata = nd, type = "prob")  
}

predictor <- Predictor$new(
  model = final_model,
  data  = X_df,
  y     = Y,
  predict.function = pred_fun_iml,
  type  = "prob"
)

set.seed(2025)
n_explain   <- min(200, nrow(X_df))
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
  summarise(mean_abs_shap = mean(abs(phi), na.rm = TRUE)) %>%
  arrange(desc(mean_abs_shap)) %>%
  slice_head(n = 10)

feature_labels <- c(
  peak_val = "Early peak value",
  zero_count = "Non-zero counts",
  c22_DN_OutlierInclude_p_001_mdrmd = "Outlier prevalence",
  c22_DN_HistogramMode_10 = "Histogram mode (bin=10)",
  c22_FC_LocalSimple_mean3_stderr = "Local fluctuation (stderr, lag=3)",
  c22_FC_LocalSimple_mean1_tauresrat = "Local autocorrelation ratio (lag=1)",
  c22_SB_BinaryStats_mean_longstretch1 = "Longest run (binary sequence)",
  c22_SB_TransitionMatrix_3ac_sumdiagcov = "Transition-matrix covariance (3-state)",
  c22_SB_MotifThree_quantile_hh = "Motif frequency (3-state, quantile)",
  c22_SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1 = "Detrended fluctuation (scale 50)"
)

top_feats <- top_feats %>%
  mutate(feature_clean = recode(feature, !!!feature_labels))

p_imp <- ggplot(top_feats, aes(x = reorder(feature_clean, mean_abs_shap), y = mean_abs_shap)) +
  geom_col(width = 0.75, fill = "#E31A1C") +
  coord_flip() +
  expand_limits(y = max(top_feats$mean_abs_shap) * 1.15) +
  labs(x = NULL, y = "Mean |SHAP|", title = "Top 10 Global Importance (mean |SHAP|)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_imp)

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
    feature = factor(feature, levels = rev(top_feats$feature)),
    feature_clean = recode(feature, !!!feature_labels)
  ) %>%
  group_by(feature, feature_clean) %>%
  mutate(
    val_pct = percent_rank(x.recoded),
    col_val = val_pct^0.6,
    zero_var = sd(x.recoded, na.rm = TRUE) == 0
  ) %>%
  ungroup()

shap_all$col_val[shap_all$zero_var] <- NA

p_bee <- ggplot(shap_all, aes(x = feature_clean, y = phi, color = col_val, alpha = col_val, size = col_val)) +
  geom_jitter(width = 0.22, height = 0, shape = 16, stroke = 0) +
  scale_color_gradientn(
    colors = c("#2166AC", "#6EA6CD", "#F7A789", "#B2182B"),
    values = c(0, 0.6, 0.9, 1),
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
