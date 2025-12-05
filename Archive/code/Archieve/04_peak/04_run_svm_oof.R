#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
proc_number <- if (length(args) >= 1) as.integer(args[1]) else 1L
if (is.na(proc_number)) proc_number <- 1L
set.seed(proc_number)

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(kernlab)   # caret::svmRadial 背后用的就是 kernlab
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

fit_and_oof_svm <- function(sp, seed = 1L) {
  set.seed(seed)
  X <- sp$X; Y <- sp$Y
  
  # 估计 sigma 的基准（用于构造合理搜索网格）
  sigma_est <- tryCatch(as.numeric(kernlab::sigest(X, frac = 1)[1]), error = function(e) NA_real_)
  if (!is.finite(sigma_est) || sigma_est <= 0) sigma_est <- 0.1
  
  ctrl <- caret::trainControl(method = "cv", number = 5, verboseIter = FALSE, allowParallel = TRUE)
  grid <- expand.grid(
    sigma = sigma_est * c(0.5, 1, 2),
    C     = c(0.5, 1, 2, 4)
  )
  
  set.seed(seed)
  svm_fit <- caret::train(
    x = X, y = Y,
    method   = "svmRadial",
    trControl= ctrl,
    tuneGrid = grid,
    metric   = "RMSE",
    preProcess = NULL
  )
  
  best_params <- svm_fit$bestTune  # 列：sigma, C
  
  folds <- caret::createFolds(Y, k = 5, list = TRUE, returnTrain = FALSE)
  
  oof_list <- lapply(seq_along(folds), function(k) {
    te_idx <- folds[[k]]
    tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
    
    # 与 caret::svmRadial 对齐：kernlab::ksvm，RBF 核（rbfdot），eps-regression
    model <- kernlab::ksvm(
      x = X[tr_idx, , drop = FALSE],
      y = Y[tr_idx],
      kernel = "rbfdot",
      kpar   = list(sigma = best_params$sigma),
      C      = best_params$C,
      type   = "eps-svr",
      scaled = FALSE
    )
    
    pred <- as.numeric(predict(model, X[te_idx, , drop = FALSE]))
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

results <- lapply(splits, fit_and_oof_svm, seed = proc_number)

dir.create("../result/svm", showWarnings = FALSE, recursive = TRUE)
outfile <- file.path("../result/svm", sprintf("svm_results_seed%03d.rds", proc_number))
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