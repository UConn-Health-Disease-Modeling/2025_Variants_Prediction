#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
proc_number <- if (length(args) >= 1) as.integer(args[1]) else 1L
if (is.na(proc_number)) proc_number <- 1L
set.seed(proc_number)

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(xgboost)
})

# ---- IO ----
feat_list    <- readRDS("feat_list.rds")
measurements <- readRDS("all_data.rds")[["measurements"]] %>%
  dplyr::select(country, lineage, days_above_30)

# ---- helper: 对齐成 X/Y/ids ----
split_one <- function(feat_df, meas_df) {
  df <- feat_df %>% dplyr::inner_join(meas_df, by = c("country","lineage"))
  id_cols   <- c("country","lineage")
  y_col     <- "days_above_30"
  feat_cols <- setdiff(names(df), c(id_cols, y_col))
  df[feat_cols] <- lapply(df[feat_cols], as.numeric)
  all_na_row <- apply(df[feat_cols], 1, function(r) all(is.na(r)))
  if (any(all_na_row)) df <- df[!all_na_row, , drop = FALSE]
  list(
    X = as.matrix(df[, feat_cols, drop = FALSE]),
    Y = as.numeric(df[[y_col]]),
    ids = df[, id_cols, drop = FALSE],
    feature_names = feat_cols
  )
}

days <- c(14, 21, 28, 35)
available <- intersect(paste0("feat_", days), names(feat_list))
splits <- setNames(
  lapply(available, function(nm) split_one(feat_list[[nm]], measurements)),
  sub("^feat_", "split_", available)
)

# ---- 5折调参 + OOF 预测 ----
fit_and_oof <- function(sp, seed = 1L) {
  X <- sp$X; Y <- sp$Y
  ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE, allowParallel = TRUE)
  grid <- expand.grid(
    nrounds = c(300, 600),
    max_depth = c(3, 5),
    eta = c(0.05, 0.10),
    gamma = c(0, 1),
    colsample_bytree = c(0.8, 1.0),
    min_child_weight = c(1, 3, 5),
    subsample = c(0.8, 1.0)
  )
  
  # grid <- expand.grid(
  #   nrounds = c(300),
  #   max_depth = c(3),
  #   eta = c(0.05),
  #   gamma = c(0),
  #   colsample_bytree = c(0.8),
  #   min_child_weight = c(1),
  #   subsample = c(0.8)
  # )
  
  set.seed(seed)
  xgb_fit <- caret::train(
    x = X, y = Y,
    method = "xgbTree",
    trControl = ctrl,
    tuneGrid  = grid,
    metric    = "RMSE",
    verbose   = FALSE
  )
  best_params <- xgb_fit$bestTune
  param <- list(
    objective        = "reg:squarederror",
    eval_metric      = "rmse",
    max_depth        = best_params$max_depth,
    eta              = best_params$eta,
    gamma            = best_params$gamma,
    colsample_bytree = best_params$colsample_bytree,
    min_child_weight = best_params$min_child_weight,
    subsample        = best_params$subsample
  )
  nrounds_best <- best_params$nrounds
  folds <- caret::createFolds(Y, k = 5, list = TRUE, returnTrain = FALSE)
  oof_list <- lapply(seq_along(folds), function(i) {
    te_idx <- folds[[i]]
    tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X[tr_idx, , drop = FALSE]), label = Y[tr_idx])
    dtest  <- xgboost::xgb.DMatrix(data = as.matrix(X[te_idx, , drop = FALSE]))
    bst <- xgboost::xgb.train(params = param, data = dtrain, nrounds = nrounds_best, verbose = 0)
    pred <- predict(bst, dtest)
    data.frame(rowIndex = te_idx, fold = i, y_true = Y[te_idx], y_pred = pred)
  })
  oof <- dplyr::bind_rows(oof_list) %>% dplyr::arrange(rowIndex)
  ids_aligned <- sp$ids[oof$rowIndex, , drop = FALSE]
  pred_df <- dplyr::bind_cols(ids_aligned, oof %>% dplyr::select(fold, y_true, y_pred))
  metrics <- caret::postResample(pred = pred_df$y_pred, obs = pred_df$y_true)
  list(
    best_params     = best_params,
    rmse            = unname(metrics["RMSE"]),
    r2              = unname(metrics["Rsquared"]),
    oof_predictions = pred_df
  )
}

results <- lapply(splits, fit_and_oof, seed = proc_number)

dir.create("../result/xgb", showWarnings = FALSE, recursive = TRUE)
outfile <- file.path("../result/xgb", sprintf("xgb_results_seed%03d.rds", proc_number))
saveRDS(results, outfile)

cat(sprintf("Saved: %s\n", outfile))
print(
  do.call(
    rbind,
    lapply(names(results), function(nm) {
      data.frame(split = nm, RMSE = results[[nm]]$rmse, R2 = results[[nm]]$r2)
    })
  )
)