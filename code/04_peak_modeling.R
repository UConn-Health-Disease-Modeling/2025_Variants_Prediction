# clear the environment
rm(list = ls())

name_map <- c(
  slope_log  = "log_slope",
  slope_y_1_7 = "short_slopes",
  longest_inc_run = "inc_run",
  zero_count = "nonzero_days",
  
  c22_DN_HistogramMode_5 = "mode_5",
  c22_DN_HistogramMode_10 = "mode_10",
  c22_CO_f1ecac = "acf_timescale",
  c22_CO_FirstMin_ac = "acf_first_min",
  c22_CO_HistogramAMI_even_2_5 = "ami2",
  c22_CO_trev_1_num = "trev",
  c22_MD_hrv_classic_pnn40 = "high_fluctuation",
  c22_SB_BinaryStats_mean_longstretch1 = "stretch_high",
  c22_SB_TransitionMatrix_3ac_sumdiagcov = "transition_variance",
  c22_PD_PeriodicityWang_th0_01 = "periodicity",
  c22_CO_Embed2_Dist_tau_d_expfit_meandiff = "embedding_dist",
  c22_IN_AutoMutualInfoStats_40_gaussian_fmmi = "ami_timescale",
  c22_FC_LocalSimple_mean1_tauresrat = "whiten_timescale",
  c22_DN_OutlierInclude_p_001_mdrmd = "outlier_timing_pos",
  c22_DN_OutlierInclude_n_001_mdrmd = "outlier_timing_neg",
  c22_SP_Summaries_welch_rect_area_5_1 = "low_freq_power",
  c22_SB_BinaryStats_diff_longstretch0 = "stretch_decreasing",
  c22_SB_MotifThree_quantile_hh = "entropy_pairs",
  c22_SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1 = "rs_range",
  c22_SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1 = "dfa",
  c22_SP_Summaries_welch_rect_centroid = "centroid_freq",
  c22_FC_LocalSimple_mean3_stderr = "forecast_error"
)

inv_logit_peak <- function(peak_share_logit, n) {
  eps <- 0.5 / (n + 1)
  peak_share_est <- plogis(peak_share_logit)  
  peak_share <- pmin(pmax(peak_share_est, eps), 1 - eps)
  return(peak_share)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)
  library(Cubist)
  library(SuperLearner)
  library(glmnet)       
  library(randomForest) 
  library(xgboost)      
  library(MASS)         
  library(gam)          
  library(kknn)         
  library(e1071)        
  library(arm)
  library(tidyr)
  library(rpart)
})

feat_list    <- readRDS("code/feat_list.rds")
measurements <- readRDS("code/all_data.rds")[["measurements"]] %>%
  dplyr::select(country, lineage, peak_share, peak_share_logit)

measurements$peak_share_cat <- cut(
  measurements$peak_share,
  breaks = c(-Inf, 0.05, 0.20, Inf),
  labels = c("<0.10", "0.10–0.20", ">0.20"),
  right = TRUE
)

# write.csv(feat_list$feat_28, "code/py_data/peak_feat_28.csv", row.names = FALSE)
# write.csv(measurements, "code/py_data/peak_measurements.csv", row.names = FALSE)




run_superlearner_binary_dual <- function(feat_data,
                                         measurements,
                                         learners = c("SL.glm", "SL.glmnet", "SL.randomForest", "SL.xgboost"),
                                         seed = 123,
                                         train_fraction = 0.5,
                                         cutoff = 0.5) {
  library(dplyr)
  library(SuperLearner)
  library(caret)
  
  set.seed(seed)
  
  # feat_data = feat_list$feat_28; learners =  c("SL.glm", "SL.glmnet", "SL.gam", "SL.rpart", "SL.randomForest", "SL.xgboost"); cutoff = .5# for test 
  
  model_df <- feat_data %>%
    dplyr::select(-auc_norm, -peak_val) %>%
    left_join(measurements, by = c("country", "lineage")) %>%
    na.omit()
  
  train_idx <- caret::createDataPartition(
    model_df$peak_share_cat,
    p = train_fraction,
    list = FALSE
  ) %>%
    as.integer()
  test_idx <- setdiff(seq_len(nrow(model_df)), train_idx)
  
  train_model_df <- model_df[train_idx, , drop = FALSE]
  test_model_df  <- model_df[test_idx, , drop = FALSE]
  
  train_df <- train_model_df %>% dplyr::select(-country, -lineage)
  test_df  <- test_model_df %>% dplyr::select(-country, -lineage)
  X_train <- train_df %>% dplyr::select(-peak_share, -peak_share_logit, -peak_share_cat)
  X_test  <- test_df %>% dplyr::select(-peak_share, -peak_share_logit, -peak_share_cat)
  
  # ---------------- Binary 1: "<0.10" vs ">0.10" ----------------
  Y1 <- ifelse(train_df$peak_share_cat == "<0.10", 0, 1)
  Y1_test <- ifelse(test_df$peak_share_cat == "<0.10", 0, 1)
  sl_fit1 <- SuperLearner(
    Y = Y1,
    X = X_train,
    SL.library = learners,
    family = binomial(),
    cvControl = list(V = 5)
  )
  pred_prob1 <- as.numeric(predict(sl_fit1, newdata = X_test, onlySL = TRUE)$pred)
  pred_class1 <- ifelse(pred_prob1 > cutoff, 1, 0)
  
  conf1 <- confusionMatrix(factor(pred_class1, levels = c(0, 1)),
                           factor(Y1_test, levels = c(0, 1)))
  acc1 <- conf1$overall["Accuracy"]
  f1_1 <- conf1$byClass["F1"]
  
  result_df1 <- test_model_df %>%
    dplyr::mutate(
      binary1_true = ifelse(Y1_test == 1, ">0.10", "<0.10"),
      binary1_pred = ifelse(pred_class1 == 1, ">0.10", "<0.10"),
      binary1_pred_prob = pred_prob1
    ) %>%
    dplyr::select(country, lineage, binary1_true, binary1_pred, binary1_pred_prob)
  
  # ---------------- Binary 2: "<0.20" vs ">0.20" ----------------
  Y2 <- ifelse(train_df$peak_share_cat == ">0.20", 1, 0)
  Y2_test <- ifelse(test_df$peak_share_cat == ">0.20", 1, 0)
  sl_fit2 <- SuperLearner(
    Y = Y2,
    X = X_train,
    SL.library = learners,
    family = binomial(),
    cvControl = list(V = 5)
  )
  pred_prob2 <- as.numeric(predict(sl_fit2, newdata = X_test, onlySL = TRUE)$pred)
  pred_class2 <- ifelse(pred_prob2 > cutoff, 1, 0)
  
  conf2 <- confusionMatrix(factor(pred_class2, levels = c(0, 1)),
                           factor(Y2_test, levels = c(0, 1)))
  acc2 <- conf2$overall["Accuracy"]
  f1_2 <- conf2$byClass["F1"]
  
  result_df2 <- test_model_df %>%
    dplyr::mutate(
      binary2_true = ifelse(Y2_test == 1, ">0.20", "<0.20"),
      binary2_pred = ifelse(pred_class2 == 1, ">0.20", "<0.20"),
      binary2_pred_prob = pred_prob2
    ) %>%
    dplyr::select(country, lineage, binary2_true, binary2_pred, binary2_pred_prob)
  
  # ---------------- Combine & reconstruct 3-class prediction ----------------
  final_df <- result_df1 %>%
    left_join(result_df2, by = c("country", "lineage")) %>%
    mutate(
      true_peak_share_cat = test_model_df$peak_share_cat,
      predicted_peak_share_cat = case_when(
        binary1_pred == "<0.10" ~ "<0.10",
        binary1_pred == ">0.10" & binary2_pred == "<0.20" ~ "0.10–0.20",
        binary1_pred == ">0.10" & binary2_pred == ">0.20" ~ ">0.20",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(country, lineage, true_peak_share_cat, predicted_peak_share_cat,
                  binary1_pred_prob, binary2_pred_prob)
  
  # ---------------- Evaluate overall 3-class performance ----------------
  final_conf <- confusionMatrix(
    factor(final_df$predicted_peak_share_cat,
           levels = c("<0.10", "0.10–0.20", ">0.20")),
    factor(final_df$true_peak_share_cat,
           levels = c("<0.10", "0.10–0.20", ">0.20"))
  )
  final_acc <- final_conf$overall["Accuracy"]
  final_f1 <- mean(final_conf$byClass[,"F1"], na.rm = TRUE)
  
  list(
    sl_fit_binary1 = sl_fit1,
    sl_fit_binary2 = sl_fit2,
    train_index = train_idx,
    test_index = test_idx,
    train_fraction = train_fraction,
    accuracy = c(binary1 = acc1, binary2 = acc2, final = final_acc),
    F1 = c(binary1 = f1_1, binary2 = f1_2, final = final_f1),
    result_df1 = result_df1,
    result_df2 = result_df2,
    final_df = final_df
  )
}



# ------------------- Setup -------------------
feature_sets <- c("feat_14", "feat_21", "feat_28", "feat_35")

learner_combos <- list(
  glm        = "SL.glm",
  bayesglm   = "SL.bayesglm",
  gam        = "SL.gam",
  elasticnet = "SL.glmnet",
  ranger     = "SL.earth",
  nnet       = "SL.nnet",
  CART       = "SL.rpart",
  ensemble   = c("SL.glm", "SL.glmnet", "SL.gam", "SL.rpart",
                 "SL.randomForest", "SL.xgboost")
)

results_list <- list()

# ------------------- Loop -------------------
for (fname in feature_sets) {
  for (mname in names(learner_combos)) {
    cat("Running SuperLearner for:", fname, "-", mname, "\n")

    res <- run_superlearner_binary_dual(
      feat_data = feat_list[[fname]],
      measurements = measurements,
      learners = learner_combos[[mname]],
      seed = 2025
    )

    results_list[[paste(fname, mname, sep = "_")]] <- res
  }
}

saveRDS(results_list, "code/peak_results_1117.rds")










results_list <- readRDS("code/peak_results_1117.rds")
summary_df <- data.frame(
  feature_set = sub("^(feat_\\d+)_.*", "\\1", names(results_list)),
  model = sub("^feat_\\d+_(.*)$", "\\1", names(results_list)),
  accuracy = sapply(results_list, function(x) {
    acc <- x$accuracy
    if (is.null(acc)) return(NA)
    if (any(grepl("final", names(acc), ignore.case = TRUE))) {
      return(as.numeric(acc[grep("final", names(acc), ignore.case = TRUE)]))
    } else {
      return(as.numeric(tail(acc, 1)))
    }
  }),
  F1 = sapply(results_list, function(x) {
    f1 <- x$F1
    if (is.null(f1)) return(NA)
    if (any(grepl("final", names(f1), ignore.case = TRUE))) {
      return(as.numeric(f1[grep("final", names(f1), ignore.case = TRUE)]))
    } else {
      return(as.numeric(tail(f1, 1)))
    }
  })
)

summary_wide <- summary_df %>%
  mutate(model = factor(model, levels = c("glm", "bayesglm", "gam", "elasticnet", "ranger", "nnet", "CART", "ensemble"))) %>%
  arrange(model, feature_set) %>%
  tidyr::pivot_wider(
    names_from = feature_set,
    values_from = c(accuracy, F1)
  ) %>%
  dplyr::select(
    model,
    accuracy_feat_14, F1_feat_14,
    accuracy_feat_21, F1_feat_21,
    accuracy_feat_28, F1_feat_28,
    accuracy_feat_35, F1_feat_35
  )

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(scales)
})

df_list <- lapply(names(results_list), function(nm) {
  if (grepl("^feat_28", nm)) {
    tmp <- results_list[[nm]]$final_df
    tmp$model <- gsub("feat_28_", "", nm)
    tmp
  } else NULL
})
df_list <- df_list[!sapply(df_list, is.null)]
all_df <- bind_rows(df_list)

country_acc <- all_df %>%
  group_by(model, country) %>%
  summarise(
    accuracy = mean(true_peak_share_cat == predicted_peak_share_cat, na.rm = TRUE),
    .groups = "drop"
  )

model_name_map <- c(
  "glm" = "GLM",
  "bayesglm" = "BayesGLM",
  "gam" = "GAM",
  "elasticnet" = "Elastic Net",
  "ranger" = "SVM",
  "nnet" = "Neural Net",
  "CART" = "CART",
  "ensemble" = "SuperLearner"
)

order_models <- c("GLM", "BayesGLM", "GAM", "Elastic Net",
                  "SVM", "Neural Net", "CART", "SuperLearner")

heat_data <- country_acc %>%
  mutate(model = sub("^feat_\\d+_", "", model)) %>%  
  mutate(model = model_name_map[model]) %>%
  mutate(model = factor(model, levels = order_models))

ggplot(
  heat_data_sim,
  aes(x = model, y = reorder(country, desc(country)), fill = accuracy_sim)
) +
  geom_tile(color = "grey90", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", accuracy_sim)),
            size = 2.8, color = "black") +
  scale_fill_gradientn(
    colors = c("#f7fbff", "#deebf7", "#9ecae1", "#3182bd", "#08519c"),
    limits = c(0.4, 1.0),
    oob = scales::squish,
    breaks = c(0.4, 0.6, 0.8, 1.0),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(x = "Model", y = "Country", fill = "Accuracy") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 11, margin = ggplot2::margin(t = 6)),
    axis.title.y = element_text(size = 11, margin = ggplot2::margin(r = 6)),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = grid::unit(4, "mm"),
    legend.key.width  = grid::unit(4, "mm"),
    legend.position = "right",
    plot.margin = ggplot2::margin(6, 6, 6, 6),
    panel.grid = element_blank()
  )

ggsave(
  filename = "PNAS/figs/peak_acc_heatmap.tiff",
  plot = last_plot(),
  width = 8.0,
  height = 5.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
