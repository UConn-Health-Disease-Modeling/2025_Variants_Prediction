#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages(library(dplyr))
source("feature_functions.R")

# Run from the code directory:
#   cd code
#   Rscript 03_extract_features.R
input_file <- "variants.rds"
train_output_file <- "train.rds"
test_output_file <- "test.rds"
reference_file <- "feat_list.rds"
input_days <- c(14L, 21L, 28L)
split_date <- as.Date("2023-05-05")

if (!requireNamespace("Rcatch22", quietly = TRUE) &&
    !requireNamespace("catch22", quietly = TRUE)) {
  stop("Install Rcatch22 before running this script to preserve the 30-column feature format.")
}

variants <- readRDS(input_file)
required_cols <- c("country", "date", "variant", "share")
missing_cols <- setdiff(required_cols, names(variants))

if (length(missing_cols) > 0L) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Enforce the workflow-wide exclusion even if an older variants.rds is supplied.
variants <- variants %>%
  filter(is.na(variant) | tolower(trimws(variant)) != "unassigned")

# Use `lineage` as the identifier so the output matches the existing feat_list format.
model_input <- variants %>%
  transmute(
    country,
    lineage = variant,
    date = as.Date(as.character(date)),
    share = as.numeric(share)
  )

# Apply only two exclusion rules: insufficient calendar span from the first
# recorded date, or a maximum share strictly above 50% in the input window.
# Missing dates inside an otherwise long-enough window represent no recorded
# sequences for that variant and are filled with share = 0.
prepare_feature_window <- function(data, input_days, dominant_threshold = 0.5) {
  data_clean <- data %>%
    group_by(country, lineage, date) %>%
    summarise(share = mean(share, na.rm = TRUE), .groups = "drop")

  periods <- data_clean %>%
    group_by(country, lineage) %>%
    summarise(
      first_date = min(date),
      last_date = max(date),
      span_days = as.integer(last_date - first_date) + 1L,
      .groups = "drop"
    )

  window_observations <- data_clean %>%
    inner_join(periods, by = c("country", "lineage")) %>%
    mutate(time = as.integer(date - first_date) + 1L) %>%
    filter(time >= 1L, time <= input_days)

  combo_flags <- periods %>%
    left_join(
      window_observations %>%
        group_by(country, lineage) %>%
        summarise(max_share_window = max(share, na.rm = TRUE), .groups = "drop"),
      by = c("country", "lineage")
    ) %>%
    mutate(
      span_long_enough = span_days >= input_days,
      exceeds_50pct = max_share_window > dominant_threshold,
      retained = span_long_enough & !exceeds_50pct
    )

  retained_periods <- combo_flags %>%
    filter(retained) %>%
    select(country, lineage, first_date)

  X <- retained_periods %>%
    tidyr::crossing(time = seq_len(input_days)) %>%
    mutate(date = first_date + time - 1L) %>%
    left_join(
      window_observations %>% select(country, lineage, date, share),
      by = c("country", "lineage", "date")
    ) %>%
    mutate(share = coalesce(share, 0)) %>%
    select(country, lineage, time, share) %>%
    tidyr::pivot_wider(
      names_from = time,
      values_from = share,
      names_prefix = "day_",
      values_fill = 0
    ) %>%
    arrange(country, lineage)

  list(X = X, combo_flags = combo_flags)
}

res_list <- setNames(
  lapply(input_days, function(d) {
    prepare_feature_window(model_input, d, dominant_threshold = 0.5)
  }),
  paste0("input_", input_days)
)

feat_list <- setNames(
  lapply(input_days, function(d) {
    extract_features_from_model_data(
      res_list[[paste0("input_", d)]]$X,
      add_auto_features = TRUE
    )
  }),
  paste0("feat_", input_days)
)

# Duration runs from the first retained date through the end of the latest run
# containing at least seven consecutive calendar days with share > 1%.
qualifying_days <- model_input %>%
  filter(share > 0.01) %>%
  distinct(country, lineage, date) %>%
  arrange(country, lineage, date) %>%
  group_by(country, lineage) %>%
  mutate(
    previous_date = lag(date),
    run_id = cumsum(is.na(previous_date) |
                      as.integer(date - previous_date) != 1L)
  ) %>%
  ungroup()

last_qualifying_runs <- qualifying_days %>%
  group_by(country, lineage, run_id) %>%
  summarise(
    run_end = max(date),
    run_length = n(),
    .groups = "drop"
  ) %>%
  filter(run_length >= 7L) %>%
  group_by(country, lineage) %>%
  summarise(last_7day_end = max(run_end), .groups = "drop")

outcomes <- model_input %>%
  group_by(country, lineage) %>%
  summarise(
    start_date = min(date),
    overall_peak_share = max(share, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  inner_join(last_qualifying_runs, by = c("country", "lineage")) %>%
  mutate(
    duration = as.integer(last_7day_end - start_date) + 1L
  )

response_list <- setNames(
  lapply(input_days, function(d) {
    initial_peaks <- res_list[[paste0("input_", d)]]$combo_flags %>%
      filter(retained) %>%
      transmute(
        country,
        lineage,
        initial_peak_share = max_share_window
      )

    feat_list[[paste0("feat_", d)]] %>%
      select(country, lineage) %>%
      inner_join(outcomes, by = c("country", "lineage")) %>%
      inner_join(initial_peaks, by = c("country", "lineage")) %>%
      transmute(
        country,
        lineage,
        duration,
        # Subsequent growth is the increase from the maximum share observed in
        # the input window to the maximum share over the complete trajectory.
        growth = 100 * (overall_peak_share - initial_peak_share)
      ) %>%
      mutate(
        duration_cat = case_when(
          duration < 30L ~ "short",
          duration <= 90L ~ "medium",
          TRUE ~ "long"
        ),
        growth_cat = case_when(
          growth <= 1 ~ "minimal",
          growth <= 5 ~ "moderate",
          TRUE ~ "large"
        ),
        duration_cat = factor(
          duration_cat,
          levels = c("short", "medium", "long"),
          ordered = TRUE
        ),
        growth_cat = factor(
          growth_cat,
          levels = c("minimal", "moderate", "large"),
          ordered = TRUE
        )
      )
  }),
  paste0("response_", input_days)
)

# Validate against the previous feature schema when the reference file exists.
if (file.exists(reference_file)) {
  reference <- readRDS(reference_file)
  reference_names <- names(reference$feat_14)
  schema_ok <- vapply(
    feat_list,
    function(x) identical(names(x), reference_names),
    logical(1)
  )
  if (!all(schema_ok)) {
    stop("Feature schema differs from ", reference_file, ".")
  }
}

response_ok <- vapply(input_days, function(d) {
  feature_keys <- feat_list[[paste0("feat_", d)]][c("country", "lineage")]
  response <- response_list[[paste0("response_", d)]]
  identical(feature_keys, response[c("country", "lineage")]) &&
    nrow(response) == nrow(feature_keys) &&
    !anyNA(response) &&
    all(response$duration >= 7L) &&
    all(response$growth >= 0) &&
    all(response$growth <= 100) &&
    identical(levels(response$duration_cat), c("short", "medium", "long")) &&
    identical(levels(response$growth_cat), c("minimal", "moderate", "large"))
}, logical(1))

if (!all(response_ok)) {
  stop("Response validation failed.")
}

# Combine predictors and categorical outcomes. Generic output names avoid
# exposing feature and outcome semantics in the shared data.
# Continuous duration and growth outcomes are intentionally not retained.
combined_data <- setNames(
  lapply(input_days, function(d) {
    features <- feat_list[[paste0("feat_", d)]]
    categorical_response <- response_list[[paste0("response_", d)]] %>%
      select(country, lineage, duration_cat, growth_cat)

    combined <- features %>%
      inner_join(categorical_response, by = c("country", "lineage")) %>%
      select(country, lineage, everything(), duration_cat, growth_cat)

    feature_count <- ncol(features) - 2L
    names(combined) <- c(
      "id_1",
      "id_2",
      sprintf("feature_%02d", seq_len(feature_count)),
      "label_1",
      "label_2"
    )

    combined
  }),
  paste0("data_", input_days)
)

# Split by each variant's first appearance across all countries. This keeps a
# lineage entirely in one dataset and prevents cross-country lineage leakage.
variant_start_dates <- model_input %>%
  group_by(lineage) %>%
  summarise(variant_start_date = min(date), .groups = "drop")

split_dataset <- function(before_split) {
  setNames(
    lapply(input_days, function(d) {
      combined_data[[paste0("data_", d)]] %>%
        left_join(variant_start_dates, by = c("id_2" = "lineage")) %>%
        filter(if (before_split) {
          variant_start_date < split_date
        } else {
          variant_start_date >= split_date
        }) %>%
        select(-variant_start_date)
    }),
    paste0(if (before_split) "train_" else "test_", input_days)
  )
}

train <- split_dataset(before_split = TRUE)
test <- split_dataset(before_split = FALSE)

validate_dataset <- function(dataset, prefix) {
  vapply(input_days, function(d) {
    split_data <- dataset[[paste0(prefix, "_", d)]]
    expected_names <- c(
      "id_1",
      "id_2",
      sprintf("feature_%02d", seq_len(ncol(feat_list[[paste0("feat_", d)]]) - 2L)),
      "label_1",
      "label_2"
    )

    identical(names(split_data), expected_names) &&
      !anyNA(split_data) &&
      is.ordered(split_data$label_1) &&
      is.ordered(split_data$label_2) &&
      identical(levels(split_data$label_1), c("short", "medium", "long")) &&
      identical(levels(split_data$label_2), c("minimal", "moderate", "large"))
  }, logical(1))
}

partition_ok <- vapply(input_days, function(d) {
  training_data <- train[[paste0("train_", d)]]
  testing_data <- test[[paste0("test_", d)]]
  all_data <- combined_data[[paste0("data_", d)]]

  nrow(training_data) + nrow(testing_data) == nrow(all_data) &&
    length(intersect(unique(training_data$id_2), unique(testing_data$id_2))) == 0L &&
    nrow(inner_join(
      training_data %>% select(id_1, id_2),
      testing_data %>% select(id_1, id_2),
      by = c("id_1", "id_2")
    )) == 0L
}, logical(1))

if (!all(validate_dataset(train, "train")) ||
    !all(validate_dataset(test, "test")) ||
    !all(partition_ok)) {
  stop("Train/test validation failed.")
}

saveRDS(train, train_output_file, compress = TRUE)
saveRDS(test, test_output_file, compress = TRUE)

for (saved_object in list(
  list(path = train_output_file, data = train),
  list(path = test_output_file, data = test)
)) {
  message(
    "Saved ", saved_object$path, ": ",
    paste(
      names(saved_object$data),
      vapply(saved_object$data, nrow, integer(1)),
      "rows x",
      vapply(saved_object$data, ncol, integer(1)),
      "columns",
      collapse = "; "
    )
  )
}
