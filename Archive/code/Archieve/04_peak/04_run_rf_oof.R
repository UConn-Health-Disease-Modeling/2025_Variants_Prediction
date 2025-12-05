#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
proc_number <- if (length(args) >= 1) as.integer(args[1]) else 1L
if (is.na(proc_number)) proc_number <- 1L
set.seed(proc_number)

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(ranger)
})

feat_list    <- readRDS("feat_list.rds")
measurements <- readRDS("all_data.rds")[["measurements"]] %>%
  dplyr::select(country, lineage, days_above_30)

split_one <- function(feat_df, meas_df) {
  df <- feat_df %>% dplyr::inner_join(meas_df, by = c("country","lineage"))
  id_cols   <- c("country","lineage")
  y_col     <- "days_above_30"
  feat_cols <- setdiff(names(df), c(id_cols, y_col))
  df[feat_cols] <- lapply(df[feat_cols], as.numeric)
  all_na_row <- apply(df[feat_cols], 1, function(r) all(is.na(r)))
  if (any(all_na_row)) df <- df[!all_na_row, , drop = FALSE]
  list(
    X = as.data.frame(df[, feat_cols, drop = FALSE]),
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

fit_and_oof_rf <- function(sp, seed = 1L) {
  X <- sp$X; Y <- sp$Y
  p <- ncol(X)
  mtry_vals <- unique(pmax(1L, round(c(sqrt(p), p/4, p/2))))
  grid <- expand.grid(
    mtry          = mtry_vals,
    splitrule     = c("variance", "extratrees"),
    min.node.size = c(1L, 5L, 10L)
  )
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    verboseIter = FALSE,
    savePredictions = "final",
    allowParallel = TRUE
  )
  set.seed(seed)
  rf_fit <- caret::train(
    x = X, y = Y,
    method = "ranger",
    trControl = ctrl,
    tuneGrid  = grid,
    metric    = "RMSE",
    preProcess = "medianImpute",
    num.trees = 1000,
    importance = "permutation",
    keep.inbag = TRUE
  )
  best_params <- rf_fit$bestTune
  pred_best <- rf_fit$pred
  for (nm in names(best_params)) {
    pred_best <- pred_best[pred_best[[nm]] == best_params[[nm]], , drop = FALSE]
  }
  oof <- pred_best[, c("rowIndex","Resample","obs","pred")]
  colnames(oof) <- c("rowIndex","fold","y_true","y_pred")
  oof <- oof[order(oof$rowIndex), ]
  ids_aligned <- sp$ids[oof$rowIndex, , drop = FALSE]
  pred_df <- dplyr::bind_cols(ids_aligned, oof[, c("fold","y_true","y_pred")])
  metrics <- caret::postResample(pred = pred_df$y_pred, obs = pred_df$y_true)
  list(
    best_params     = best_params,
    rmse            = unname(metrics["RMSE"]),
    r2              = unname(metrics["Rsquared"]),
    oof_predictions = pred_df
  )
}

results <- lapply(splits, fit_and_oof_rf, seed = proc_number)

dir.create("../result/rf", showWarnings = FALSE, recursive = TRUE)
outfile <- file.path("../result/rf", sprintf("rf_results_seed%03d.rds", proc_number))
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