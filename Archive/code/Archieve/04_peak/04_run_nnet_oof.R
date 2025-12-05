#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
proc_number <- if (length(args) >= 1) as.integer(args[1]) else 1L
if (is.na(proc_number)) proc_number <- 1L
set.seed(proc_number)

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(nnet)
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

fit_and_oof_nnet <- function(sp, seed = 1L) {
  X <- sp$X; Y <- sp$Y
  set.seed(seed)
  
  ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE, allowParallel = TRUE)
  
  grid <- expand.grid(
    size  = c(5, 10, 20),     # 隐层神经元数
    decay = c(0, 1e-4, 1e-3, 1e-2)  # L2 正则
  )
  
  # 这里使用 medianImpute + center + scale，均在每个 CV 训练折内拟合（避免信息泄漏）
  nnet_fit <- caret::train(
    x = X, y = Y,
    method = "nnet",
    trControl = ctrl,
    tuneGrid  = grid,
    metric    = "RMSE",
    preProcess = c("medianImpute", "center", "scale"),
    linout = TRUE,    # 回归
    trace  = FALSE,
    maxit  = 500,
    MaxNWts = max(10000, 10 * ncol(X) + 10)
  )
  
  best_params <- nnet_fit$bestTune
  
  folds <- caret::createFolds(Y, k = 5, list = TRUE, returnTrain = FALSE)
  oof_list <- lapply(seq_along(folds), function(k) {
    te_idx <- folds[[k]]
    tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
    
    # 按最佳超参训练该折模型；preProcess 同样在训练折内拟合并作用于测试折
    m <- caret::train(
      x = X[tr_idx, , drop = FALSE],
      y = Y[tr_idx],
      method = "nnet",
      trControl = trainControl(method = "none"),
      tuneGrid  = best_params,
      preProcess = c("medianImpute", "center", "scale"),
      linout = TRUE,
      trace  = FALSE,
      maxit  = 500,
      MaxNWts = max(10000, 10 * ncol(X) + 10)
    )
    
    pred <- as.numeric(predict(m, newdata = X[te_idx, , drop = FALSE]))
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

results <- lapply(splits, fit_and_oof_nnet, seed = proc_number)

dir.create("../result/nnet", showWarnings = FALSE, recursive = TRUE)
outfile <- file.path("../result/nnet", sprintf("nnet_results_seed%03d.rds", proc_number))
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