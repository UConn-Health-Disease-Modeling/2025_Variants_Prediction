#!/usr/bin/env Rscript
# setwd("code/")
# setwd(here::here()) 

args <- commandArgs(trailingOnly = TRUE)
proc_number <- if (length(args) >= 1) as.integer(args[1]) else 1L
if (is.na(proc_number)) proc_number <- 1L
set.seed(proc_number)

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(Cubist)
  library(ranger)
})

feat_list    <- readRDS("../feat_list.rds")
measurements <- readRDS("../all_data.rds")[["measurements"]] %>%
  dplyr::select(country, lineage, days_above_10_cat)

split_one <- function(feat_df, meas_df,
                      y_col = "days_above_10_cat",
                      id_cols = c("country","lineage"),
                      ordered_levels = c("0","1–30","31–100","100+"),
                      make_ordinal = FALSE) {
  df <- dplyr::inner_join(feat_df, meas_df, by = id_cols)
  feat_cols <- setdiff(names(df), c(id_cols, y_col))
  df[feat_cols] <- lapply(df[feat_cols], function(x) {
    if (is.numeric(x)) x else suppressWarnings(as.numeric(x))
  })
  all_na_row <- apply(df[feat_cols], 1, function(r) all(is.na(r)))
  if (any(all_na_row)) df <- df[!all_na_row, , drop = FALSE]
  y_raw <- df[[y_col]]
  if (is.numeric(y_raw)) y_raw <- as.character(y_raw)
  observed <- unique(as.character(y_raw))
  use_levels <- if (length(intersect(ordered_levels, observed)) >= 2) {
    ordered_levels[ordered_levels %in% observed]
  } else {
    sort(observed)
  }
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

days <- c(14, 21, 28, 35)
available <- intersect(paste0("feat_", days), names(feat_list))
splits <- setNames(
  lapply(available, function(nm) {
    split_one(
      feat_df = feat_list[[nm]],
      meas_df = measurements,
      y_col = "days_above_10_cat",
      ordered_levels = c("0","1–30","31–100","100+"),
      make_ordinal = TRUE
    )
  }),
  sub("^feat_", "split_", available)
)

splits1_countries <- c(
  "United Kingdom","United States",
  "Germany","Japan","Canada","France"
)

splits2_countries <- c(
  "Denmark", "Sweden","India","Brazil",
  "Australia","Spain","Italy","Austria","South Korea"
)
make_country_splits <- function(splits, countries) {
  lapply(splits, function(sp) {
    keep <- sp$ids$country %in% countries
    lapply(sp, function(x) {
      if (is.data.frame(x) || is.matrix(x)) x[keep, , drop = FALSE] else x[keep]
    })
  })
}

splits1 <- make_country_splits(splits, splits1_countries)
splits2 <- make_country_splits(splits, splits2_countries)

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

fit_and_oof_classif <- function(sp, seed = 1L) {
  X <- sp$X
  Y <- if (!is.null(sp$Y_factor)) sp$Y_factor else as.factor(sp$Y)
  if (!is.factor(Y)) Y <- as.factor(Y)
  if (nlevels(Y) < 2) stop("Y has fewer than 2 classes.")
  set.seed(seed)
  
  keep_global <- colSums(!is.na(X)) > 0
  Xg <- X[, keep_global, drop = FALSE]
  
  ctrl <- caret::trainControl(
    method = "cv",
    number = 5,
    classProbs = FALSE,
    summaryFunction = defaultSummary,
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  
  grid <- expand.grid(
    mtry = unique(pmax(1, round(ncol(Xg) * c(0.1, 0.25, 0.5)))),
    splitrule = "gini",
    min.node.size = c(1, 5, 10)
  )
  
  ranger_fit <- caret::train(
    x = as.data.frame(Xg),
    y = Y,
    method = "ranger",
    trControl = ctrl,
    tuneGrid  = grid,
    metric    = "Accuracy",
    num.trees = 500,
    importance = "impurity",
    preProcess = c("medianImpute","zv")
  )
  
  best_params <- ranger_fit$bestTune
  
  folds <- caret::createFolds(Y, k = 5, list = TRUE, returnTrain = FALSE)
  
  oof_list <- lapply(seq_along(folds), function(k) {
    te_idx <- folds[[k]]
    tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
    
    Xi_tr_raw <- X[tr_idx, , drop = FALSE]
    Xi_te_raw <- X[te_idx, , drop = FALSE]
    
    non_all_na <- colSums(!is.na(Xi_tr_raw)) > 0
    Xi_tr_raw <- Xi_tr_raw[, non_all_na, drop = FALSE]
    Xi_te_raw <- Xi_te_raw[, non_all_na, drop = FALSE]
    
    imp_tr <- impute_median(Xi_tr_raw, med = NULL)
    keep_med <- is.finite(imp_tr$med)
    Xtr_i <- imp_tr$X[, keep_med, drop = FALSE]
    med_used <- imp_tr$med[keep_med]
    
    imp_te <- impute_median(Xi_te_raw[, keep_med, drop = FALSE], med = med_used)
    Xte_i <- imp_te$X
    
    nzv_idx <- caret::nearZeroVar(Xtr_i, saveMetrics = FALSE)
    if (length(nzv_idx)) {
      Xtr_i <- Xtr_i[, -nzv_idx, drop = FALSE]
      Xte_i <- Xte_i[, -nzv_idx, drop = FALSE]
    }
    
    mdl <- ranger::ranger(
      x = as.data.frame(Xtr_i),
      y = Y[tr_idx],
      mtry = best_params$mtry,
      splitrule = as.character(best_params$splitrule),
      min.node.size = best_params$min.node.size,
      num.trees = 500,
      respect.unordered.factors = "order"
    )
    
    pr <- predict(mdl, data = as.data.frame(Xte_i))$predictions
    if (is.matrix(pr)) pr <- colnames(pr)[max.col(pr, ties.method = "first")]
    
    data.frame(rowIndex = te_idx, fold = k, y_true = Y[te_idx], y_pred = factor(pr, levels = levels(Y)))
  })
  
  oof <- dplyr::bind_rows(oof_list) %>% dplyr::arrange(rowIndex)
  ids_aligned <- sp$ids[oof$rowIndex, , drop = FALSE]
  pred_df <- dplyr::bind_cols(ids_aligned, oof %>% dplyr::select(fold, y_true, y_pred))
  mets <- macro_metrics(truth = pred_df$y_true, pred = pred_df$y_pred)
  
  list(
    best_params       = best_params,
    accuracy          = unname(mets["accuracy"]),
    macro_f1          = unname(mets["macro_f1"]),
    balanced_accuracy = unname(mets["balanced_accuracy"]),
    oof_predictions   = pred_df
  )
}

results1 <- lapply(splits1, fit_and_oof_classif, seed = proc_number)
results2 <- lapply(splits2, fit_and_oof_classif, seed = proc_number)

results <- list(
  results1 = results1, 
  results2 = results2
)

dir.create("../../result/Share_30/cubist", showWarnings = FALSE, recursive = TRUE)
outfile <- file.path("../../result/Share_30/cubist", sprintf("cubist_results_seed%03d.rds", proc_number))
saveRDS(results, outfile)