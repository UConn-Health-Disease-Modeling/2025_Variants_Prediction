#!/usr/bin/env Rscript

# Print country-level lineage diversity and total sequence counts from the GISAID
# daily lineage summary. "Unassigned" lineages are excluded before any country
# counts or rankings are calculated. A lineage is counted when it has
# numerator > 0, and sum_cases is the sum of those positive numerator values.

suppressPackageStartupMessages(library(dplyr))

raw_data_relative_path <- file.path(
  "reference",
  "UKHSA-UConn-variant-modelling",
  "variant_modelling",
  "data",
  "summary_GISAID_20240918.csv"
)

find_project_root <- function() {
  search_starts <- getwd()

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- sub("^--file=", "", file_arg[1L])
    search_starts <- c(dirname(normalizePath(script_path)), search_starts)
  }

  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    editor_path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )
    if (nzchar(editor_path)) {
      search_starts <- c(dirname(normalizePath(editor_path)), search_starts)
    }
  }

  for (start in unique(search_starts)) {
    current <- normalizePath(start)
    repeat {
      if (file.exists(file.path(current, raw_data_relative_path))) {
        return(current)
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }

  stop("Could not locate the project root containing: ", raw_data_relative_path)
}

project_root <- find_project_root()
raw_data_file <- file.path(project_root, raw_data_relative_path)
raw_data <- read.csv(raw_data_file, check.names = FALSE)

required_columns <- c("country", "lineage", "numerator")
missing_columns <- setdiff(required_columns, names(raw_data))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

country_summary <- raw_data %>%
  filter(
    !is.na(country), nzchar(trimws(country)),
    !is.na(lineage), nzchar(trimws(lineage)),
    tolower(trimws(lineage)) != "unassigned",
    !is.na(numerator), numerator > 0
  ) %>%
  mutate(country = recode(country, "USA" = "United States")) %>%
  group_by(country) %>%
  summarise(
    n_lineages = n_distinct(lineage),
    total_sequences = sum(numerator),
    .groups = "drop"
  ) %>%
  arrange(desc(n_lineages), desc(total_sequences), country) %>%
  mutate(rank = row_number()) %>%
  select(rank, country, n_lineages, total_sequences)

print(country_summary, n = Inf)

top15_countries <- country_summary %>%
  slice_head(n = 15L)

top15_output_file <- file.path(getwd(), "top15_countries.csv")
write.csv(top15_countries, top15_output_file, row.names = FALSE)

message("Saved Top 15 country list: ", normalizePath(top15_output_file))
