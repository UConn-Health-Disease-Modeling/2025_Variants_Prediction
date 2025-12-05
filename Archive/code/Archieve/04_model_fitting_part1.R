# ------------------------------------------------------------------------------
# 04_model_fitting_part1.R
# Usage:
#   Rscript 04_model_fitting_part1.R 1
#   Rscript 04_model_fitting_part1.R 7   # seed=7
# Output:
#   result/tune/tune_proc<proc_id>.rds
# ------------------------------------------------------------------------------

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(caret)
  library(glmnet)
})

# your modeling helpers (must provide fit_eval_models)
source("code/functions_ml_utils.R")

# ------------------------------- I/O ------------------------------------------
all_data      <- readRDS("code/all_data.rds")
features_list <- readRDS("code/features_list.rds")   # list: day14/day21/day28/day35
vol_groups    <- readRDS("code/vol_groups_view.rds") # columns: country, volatility_group (factor/chr)
measurements  <- all_data$measurements               # columns: country, classified_label, auc_raw, peak_share, days_above_30

out_dir <- "result/tune"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------- proc_id -> seed & out ------------------------------
args <- commandArgs(trailingOnly = TRUE)
proc_id <- suppressWarnings(as.integer(args[1]))
if (is.na(proc_id) || proc_id <= 0) {
  stop("Please provide a positive integer proc_id, e.g. `Rscript 04_model_fitting_part1.R 1`.")
}
seed <- proc_id
outfile <- file.path(out_dir, sprintf("tune_proc%d.rds", proc_id))

message(sprintf("==> proc_id = %d  (seed = %d)", proc_id, seed))
message("Results will be saved to: ", outfile)

# --------------------------- helpers ------------------------------------------
# 构造每个 day 的按组对齐的 X/Y
build_pairs_by_group <- function(X_all, Y_all, vol_groups) {
  X_all <- X_all %>% dplyr::left_join(vol_groups, by = "country")
  
  Y_all <- Y_all %>%
    dplyr::left_join(vol_groups, by = "country") %>%
    # keep only keys present in X
    dplyr::filter(paste(country, classified_label) %in%
                    paste(X_all$country, X_all$classified_label)) %>%
    # align order to X
    dplyr::slice(match(paste(X_all$country, X_all$classified_label),
                       paste(country, classified_label)))
  
  groups <- X_all %>%
    dplyr::pull(volatility_group) %>%
    as.character() %>%
    unique() %>%
    sort()
  
  purrr::map(setNames(groups, groups), function(g) {
    Xg <- X_all %>% dplyr::filter(as.character(volatility_group) == g)
    Yg <- Y_all %>%
      dplyr::filter(as.character(volatility_group) == g) %>%
      dplyr::filter(paste(country, classified_label) %in%
                      paste(Xg$country, Xg$classified_label)) %>%
      dplyr::slice(match(paste(Xg$country, Xg$classified_label),
                         paste(country, classified_label)))
    list(X = Xg, Y = Yg)
  })
}

# --------------------------- main loop ----------------------------------------
day_keys <- c("day14","day21","day28","day35")
targets  <- c("peak_share", "auc_raw", "days_above_30")

all_train_metrics <- list()
all_test_metrics  <- list()
rid <- 0L

for (day_key in day_keys) {
  message("---- Day set: ", day_key, " ----")
  
  X_all <- features_list[[day_key]]
  if (is.null(X_all)) {
    warning("features_list[[", day_key, "]] is NULL. Skipping.")
    next
  }
  Y_all <- measurements %>%
    dplyr::select(country, classified_label, auc_raw, peak_share, days_above_30)
  
  pairs_by_group <- build_pairs_by_group(X_all, Y_all, vol_groups)
  groups <- intersect(c("1","2","3","4"), names(pairs_by_group))
  
  for (g in groups) {
    Xg <- pairs_by_group[[g]]$X
    Yg <- pairs_by_group[[g]]$Y
    
    for (tgt in targets) {
      message(sprintf(">> [%s] group=%s, target=%s, seed=%d", day_key, g, tgt, seed))
      
      if (!tgt %in% names(Yg)) {
        warning(sprintf("Target '%s' not in Y (group %s). Skipping.", tgt, g))
        next
      }
      y_vec <- as.numeric(Yg[[tgt]])
      
      # 仅抑制警告（如 xgboost 的 ntree_limit），不屏蔽 error
      fit <- tryCatch(
        suppressWarnings( fit_eval_models(Xg, y_vec, seed = seed) ),
        error = function(e) {
          warning(sprintf("[day=%s | group=%s | target=%s] %s", day_key, g, tgt, e$message))
          NULL
        }
      )
      if (is.null(fit)) next
      
      rid <- rid + 1L
      all_train_metrics[[rid]] <- fit$train_metrics %>%
        dplyr::mutate(day = day_key, group = g, target = tgt, seed = seed, .before = 1)
      all_test_metrics[[rid]] <- fit$test_metrics %>%
        dplyr::mutate(day = day_key, group = g, target = tgt, seed = seed, .before = 1)
    }
  }
}

train_metrics_df <- dplyr::bind_rows(all_train_metrics)
test_metrics_df  <- dplyr::bind_rows(all_test_metrics)

saveRDS(
  list(
    proc_id         = proc_id,
    seed            = seed,
    train_metrics   = train_metrics_df,
    test_metrics    = test_metrics_df
  ),
  file = outfile
)

message("Saved: ", outfile)
message("Train rows: ", nrow(train_metrics_df), " | Test rows: ", nrow(test_metrics_df))