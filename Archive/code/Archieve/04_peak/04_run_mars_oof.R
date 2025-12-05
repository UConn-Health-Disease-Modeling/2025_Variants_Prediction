#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
proc_number <- if (length(args) >= 1) as.integer(args[1]) else 1L
if (is.na(proc_number)) proc_number <- 1L
set.seed(proc_number)

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(earth)   # MARS
})

# -------------------------
# 读取数据
# -------------------------
feat_list    <- readRDS("feat_list.rds")
measurements <- readRDS("all_data.rds")[["measurements"]] %>%
  dplyr::select(country, lineage, days_above_30)

# -------------------------
# 对齐到 X / Y / ids
# -------------------------
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

# -------------------------
# 简单中位数填充（返回填充后的矩阵及使用的中位数向量）
# -------------------------
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

# -------------------------
# 5 折调参 + OOF 预测（MARS）
# -------------------------
fit_and_oof_mars <- function(sp, seed = 1L) {
  X <- sp$X; Y <- sp$Y
  set.seed(seed)
  
  # caret 调参：earth 的关键参数是 degree（交互阶数）和 nprune（剪枝后基函数数）
  ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE, allowParallel = TRUE)
  
  # 先对全体做一次中位数填充，仅用于 caret::train 的一致性（避免 NA）
  imp_all <- impute_median(X, med = NULL)
  X_all_i <- imp_all$X
  
  # 为了给 nprune 提供合适范围，先训练一个粗模型估计最大可用基函数数目
  tmp_fit <- earth::earth(x = X_all_i, y = Y, degree = 2, trace = 0)
  max_terms <- length(tmp_fit$selected.terms)
  # nprune 的搜索网格不要超过可选项
  npr_grid <- sort(unique(pmin(c(10, 20, 30, 40, 50, 75, 100), max_terms)))
  if (length(npr_grid) == 0) npr_grid <- min(10L, max_terms)
  
  grid <- expand.grid(
    degree = c(1, 2),      # 只到二阶交互，稳妥
    nprune = npr_grid
  )
  
  mars_fit <- caret::train(
    x = X_all_i, y = Y,
    method = "earth",
    trControl = ctrl,
    tuneGrid  = grid,
    metric    = "RMSE",
    # 其余参数可传入 ... 到 earth()：例如 penalty, nk 等
    # 这里保持默认以减少超参维度
  )
  
  best_params <- mars_fit$bestTune   # 列：degree, nprune
  
  # 基于最佳参数做 5 折 OOF 预测（每折训练时仅用训练折的中位数做填充，避免信息泄漏）
  folds <- caret::createFolds(Y, k = 5, list = TRUE, returnTrain = FALSE)
  oof_list <- lapply(seq_along(folds), function(k) {
    te_idx <- folds[[k]]
    tr_idx <- setdiff(seq_len(nrow(X)), te_idx)
    
    imp_tr <- impute_median(X[tr_idx, , drop = FALSE], med = NULL)
    Xtr_i  <- imp_tr$X
    imp_te <- impute_median(X[te_idx, , drop = FALSE],  med = imp_tr$med)
    Xte_i  <- imp_te$X
    
    m <- earth::earth(
      x = Xtr_i, y = Y[tr_idx],
      degree = best_params$degree,
      nprune = best_params$nprune,
      trace  = 0
    )
    pred <- as.numeric(predict(m, newdata = Xte_i))
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

results <- lapply(splits, fit_and_oof_mars, seed = proc_number)

dir.create("../result/mars", showWarnings = FALSE, recursive = TRUE)
outfile <- file.path("../result/mars", sprintf("mars_results_seed%03d.rds", proc_number))
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