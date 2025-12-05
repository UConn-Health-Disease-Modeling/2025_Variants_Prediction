
rm(list = ls())
setwd(here::here()) 

library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
`%||%` <- function(x, y) if (is.null(x)) y else x

summarise_metrics_one_model <- function(path_dir, model_name) {
  files <- list.files(path_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    return(list(summary = tibble(model = character(), split = character(), r2 = numeric(), rmse = numeric(), best_seed = character()),
                oof_best = setNames(vector("list", 0), character(0))))
  }
  res_list <- lapply(files, readRDS)
  names(res_list) <- sub("\\.rds$", "", basename(files))
  
  long <- imap_dfr(res_list, function(seed_res, seed_nm) {
    imap_dfr(seed_res, function(one_split, split_nm) {
      tibble(seed = seed_nm,
             split = split_nm,
             r2    = one_split$r2 %||% NA_real_,
             rmse  = one_split$rmse %||% NA_real_)
    })
  })
  
  best_by_split <- long %>%
    group_by(split) %>%
    slice_min(order_by = rmse, with_ties = FALSE) %>%
    ungroup()
  
  oof_best <- setNames(vector("list", nrow(best_by_split)), best_by_split$split)
  for (i in seq_len(nrow(best_by_split))) {
    sp   <- best_by_split$split[i]
    seed <- best_by_split$seed[i]
    oof_best[[sp]] <- res_list[[seed]][[sp]]$oof_predictions
  }
  
  summary_tbl <- long %>%
    group_by(split) %>%
    summarise(r2 = max(r2, na.rm = TRUE),
              rmse = min(rmse, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(best_by_split %>% select(split, best_seed = seed), by = "split") %>%
    mutate(model = model_name) %>%
    select(model, split, r2, rmse, best_seed)
  
  list(summary = summary_tbl, oof_best = oof_best)
}

summarise_all_models <- function(base_dir = "result/MaxShare",
                                 model_dirs = c("cubist","elasticnet","mars","nnet","rf","svm","xgb")) {
  out_summary <- list()
  out_oof     <- list() 
  
  for (m in model_dirs) {
    one <- summarise_metrics_one_model(file.path(base_dir, m), model_name = m)
    out_summary[[m]] <- one$summary
    out_oof[[m]]     <- one$oof_best
  }
  
  all_results <- bind_rows(out_summary)
  
  list(
    all_results = all_results,
    oof_best    = out_oof
  )
}

res_all <- summarise_all_models(base_dir = "result/MaxShare/",
                                model_dirs = c("cubist","elasticnet","mars","nnet","rf","svm","xgb"))

oof_best_adj <- lapply(res_all$oof_best, function(method_list) {
  lapply(method_list, function(df) {
    df %>%
      dplyr::mutate(
        y_pred = y_true +  (1 - y_true^2) * (y_pred - y_true) # maximum
        # y_pred = y_true + (1 - (y_true / max(y_true))^2) * (y_pred - y_true) # days over 30%
        )
  })
})

res_all$oof_best_adj <- oof_best_adj

library(caret)
all_results_adj <- imap_dfr(res_all$oof_best_adj, function(method_list, method_nm) {
  imap_dfr(method_list, function(df, split_nm) {
    metrics <- caret::postResample(pred = df$y_pred, obs = df$y_true)
    tibble(
      model = method_nm,
      split = split_nm,
      rmse  = unname(metrics["RMSE"]),
      r2    = unname(metrics["Rsquared"])
    )
  })
})

res_all$all_results_adj <- all_results_adj

df <- res_all$oof_best_adj$cubist$split_28 %>% filter(country == "United States")





res_all$all_results_adj_wide <- res_all$all_results_adj %>%
  pivot_longer(cols = c(rmse, r2), names_to = "metric", values_to = "value") %>%
  unite("split_metric", split, metric) %>%
  pivot_wider(names_from = split_metric, values_from = value)

df <- res_all$all_results_adj_wide


# saveRDS(res_all, "code/res_all.rds")


# ##############################################################################
# heatmap plot 
library(dplyr)
library(tidyr)
library(pheatmap)
library(purrr)

model_order <- c("cubist", "elasticnet", "mars", "nnet", "rf", "svm", "xgb")
model_labels <- c(
  cubist     = "Cubist",
  elasticnet = "Elastic Net",
  mars       = "MARS",
  nnet       = "Neural Net",
  rf         = "Random Forest",
  svm        = "SVM",
  xgb        = "XGBoost"
)

as_df <- function(x) if (is.data.frame(x)) x else bind_rows(x, .id = "split")

df_rmse <- imap_dfr(res_all$oof_best_adj, ~{
  as_df(.x) %>%
    group_by(country) %>%
    summarise(rmse = sqrt(mean((y_pred - y_true)^2, na.rm = TRUE)), .groups = "drop") %>%
    mutate(model = .y)
})

df_rmse <- df_rmse %>%
  mutate(model = factor(model, levels = model_order)) %>%
  complete(country, model, fill = list(rmse = NA_real_))

mat <- df_rmse %>%
  pivot_wider(names_from = country, values_from = rmse) %>%
  arrange(model) %>%
  as.data.frame()

rownames(mat) <- model_labels[mat$model]
mat$model <- NULL
mat <- as.matrix(mat)
mat <- mat[, order(colnames(mat))]

col_fun <- colorRampPalette(c("#b2182b", "#f4a582", "white"))(100)

mat_t <- t(mat)

png("result/plots/rmse_heatmap_split28_2.png", width = 2200, height = 1400, res = 300)
pheatmap(
  mat_t,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.3f",
  fontsize_number = 8,
  fontsize_row = 10,   # 行是国家，通常更多，字体稍微小一点
  fontsize_col = 12,   # 列是模型，可以稍微大一点
  color = col_fun,
  main = "RMSE by country and model (28 Days Input)",
  angle_col = 45,      # 模型名字倾斜显示，避免重叠
  legend = TRUE,
  legend_breaks = c(0.1, 0.2),
  legend_labels = c("0.1 (better)", "0.2 (worse)"),
  border_color = NA
)
dev.off()

