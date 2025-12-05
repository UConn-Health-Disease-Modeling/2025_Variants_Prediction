#!/usr/bin/env Rscript
# setwd("code/04_duration/") 

args <- commandArgs(trailingOnly = TRUE)
proc_number <- if (length(args) >= 1) as.integer(args[1]) else 1L
if (is.na(proc_number)) proc_number <- 1L
set.seed(proc_number)

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(glmnet)
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

days <- c(14, 21, 28, 35)
available <- intersect(paste0("feat_", days), names(feat_list))
splits <- setNames(
  lapply(available, function(nm) split_one(feat_list[[nm]], measurements)),
  sub("^feat_", "split_", available)
)

splits1_countries <- c(
  "United Kingdom","United States",
  "Germany","Japan","Canada","France"
)
splits2_countries <- c(
  "Denmark","Sweden","India","Brazil",
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
  oof <- oof[oof$alpha == best_params$alpha & abs(oof$lambda - best_params$lambda) < .Machine$double.eps^0.5, , drop = FALSE]
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
    varimp            = caret::varImp(en_fit)
  )
}

results1 <- lapply(splits1, fit_and_oof_elasticnet, seed = proc_number)
results2 <- lapply(splits2, fit_and_oof_elasticnet, seed = proc_number)

results <- list(results1 = results1, results2 = results2)

dir.create("../../result/Share_30/elasticnet", showWarnings = FALSE, recursive = TRUE)
outfile <- file.path("../../result/Share_30/elasticnet", sprintf("elasticnet_results_seed%03d.rds", proc_number))
saveRDS(results, outfile)