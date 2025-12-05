rm(list = ls())

library(dplyr)
library(purrr)
library(readr)

# set folder where results are saved
res_dir <- "result/tune"

# proc numbers you ran
procs <- 1:5

# read all RDS files
all_results <- map(procs, function(p) {
  f <- file.path(res_dir, sprintf("tune_proc%d.rds", p))
  if (!file.exists(f)) {
    warning(sprintf("File not found: %s", f))
    return(NULL)
  }
  out <- readRDS(f)
  # make sure each has test/train metrics
  list(
    proc         = p,
    train_metrics = out$train_metrics,
    test_metrics  = out$test_metrics
  )
})

# drop any NULLs
all_results <- compact(all_results)

# combine into data frames with proc id
train_df <- bind_rows(lapply(all_results, function(x) mutate(x$train_metrics, proc = x$proc)))
test_df  <- bind_rows(lapply(all_results, function(x) mutate(x$test_metrics,  proc = x$proc)))


best_results <- test_df %>%
  dplyr::group_by(day, group, target, model) %>%
  dplyr::arrange(dplyr::desc(tidyr::replace_na(R2, -Inf)), RMSE, .by_group = TRUE) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(day, group, target, model) %>%
  dplyr::select(-seed, -proc) %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::where(is.numeric),
      ~ round(.x, 2)
    )
  )

best_results <- best_results %>%
  dplyr::mutate(
    RMSE = RMSE / 2,
    MSE  = MSE / 4
  )

