#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
proc_number <- if (length(args) >= 1) as.integer(args[1]) else 1L
if (is.na(proc_number)) proc_number <- 1L
set.seed(proc_number)

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(Cubist)
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

fit_and_oof_cubist <- function(sp, seed = 1L) {
  X <- sp$X; Y <- sp$Y
  set.seed(seed)
  ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE, allowParallel = TRUE)
  grid <- expand.grid(committees = c(10, 50, 100),
                      neighbors  = c(0, 3, 5))
  xgb_fit <- caret::train(
    x = X, y = Y,
    method = "cubist",
    trControl = ctrl,
    tuneGrid  = grid,
    metric    = "RMSE"
  )
  best_params <- xgb_fit$bestTune
  folds <- caret::createFolds(Y, k = 5, list = TRUE, returnTrain = FALSE)
  
  oof_list <- lapply(seq_along(folds), function(k) {
    te_idx <- folds[[k]]
    tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
    
    imp_tr <- impute_median(X[tr_idx, , drop = FALSE], med = NULL)
    Xtr_i  <- imp_tr$X
    imp_te <- impute_median(X[te_idx, , drop = FALSE], med = imp_tr$med)
    Xte_i  <- imp_te$X
    
    mdl <- Cubist::cubist(
      x = as.data.frame(Xtr_i),
      y = Y[tr_idx],
      committees = best_params$committees,
      neighbors  = best_params$neighbors
    )
    pred <- predict(mdl, newdata = as.data.frame(Xte_i))
    data.frame(rowIndex = te_idx, fold = k, y_true = Y[te_idx], y_pred = pred)
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

results <- lapply(splits, fit_and_oof_cubist, seed = proc_number)

dir.create("../result/cubist", showWarnings = FALSE, recursive = TRUE)
outfile <- file.path("../result/cubist", sprintf("cubist_results_seed%03d.rds", proc_number))
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