#!/usr/bin/env Rscript

# Prepare comparable variant trajectories for the 15 selected countries.
#
# Workflow:
#   1. Read every available date for the countries in top15_countries.csv and
#      exclude records whose lineage is "Unassigned".
#   2. Expand reversible Pango aliases and harmonize reporting resolution across
#      countries. This changes labels only; sequences are not reassigned.
#   3. Add a historical WHO label when one can be resolved. Unmatched variants
#      remain in the data with who_variant = NA.
#   4. Apply one eligibility rule to each country-variant trajectory:
#        peak share > 1% AND at least 7 consecutive calendar days with share > 1%.
#   5. Start each eligible trajectory on the first day of its earliest
#      qualifying seven-day run.
#   6. Save the trimmed trajectories and diagnostics, then print a compact count
#      table for every main processing step.

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

raw_data_relative_path <- file.path(
  "reference",
  "UKHSA-UConn-variant-modelling",
  "variant_modelling",
  "data",
  "summary_GISAID_20240918.csv"
)

# Run from the code directory:
#   cd code
#   Rscript 02_select_variants.R
project_root <- normalizePath("..", mustWork = TRUE)
raw_data_file <- file.path(project_root, raw_data_relative_path)
country_file <- file.path(project_root, "code", "top15_countries.csv")
output_dir <- file.path(project_root, "result", "02_group_variants")
output_rds_file <- file.path(project_root, "code", "variants.rds")

expected_pangolin_version <- "4.4"
expected_pangolin_data_version <- "1.39"
minimum_nonrecombinant_depth <- 2L  # permits B.1.1, excludes B and B.1
share_threshold <- 0.01
minimum_run_days <- 7L
diagnostic_high_share_threshold <- 0.20

who_reference <- data.table(
  who_variant = c(
    "Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Epsilon", "Zeta",
    "Eta", "Theta", "Iota", "Kappa", "Lambda", "Mu", "Omicron"
  ),
  reported_root = c(
    "B.1.1.7", "B.1.351", "P.1", "B.1.617.2", "B.1.427", "B.1.429",
    "P.2", "B.1.525", "P.3", "B.1.526", "B.1.617.1", "C.37",
    "B.1.621", "B.1.1.529"
  ),
  canonical_root = c(
    "B.1.1.7", "B.1.351", "B.1.1.28.1", "B.1.617.2", "B.1.427",
    "B.1.429", "B.1.1.28.2", "B.1.525", "B.1.1.28.3", "B.1.526",
    "B.1.617.1", "B.1.1.1.37", "B.1.621", "B.1.1.529"
  )
)

pangolin_bin <- Sys.getenv(
  "PANGOLIN_BIN",
  unname(Sys.which("pangolin"))
)
# Pangolin calls UShER, gofasta, and minimap2 as subprocesses, so the frozen
# environment's bin directory must also be visible when the R process itself
# was not launched through `micromamba run`.
if (nzchar(pangolin_bin)) {
  Sys.setenv(PATH = paste(
    dirname(pangolin_bin),
    Sys.getenv("PATH"),
    sep = .Platform$path.sep
  ))
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

fail <- function(...) stop(paste0(...), call. = FALSE)

# Return the three counts used to monitor data retention through the workflow.
# `cases_sequences` is the sum of numerator, i.e. submitted sequences rather
# than epidemiological case counts.
summarize_step <- function(step, data, variant_col) {
  if (!all(c("country", variant_col, "numerator") %in% names(data))) {
    fail("Cannot summarize step '", step, "': required columns are missing.")
  }

  country_variant_pairs <- unique(data[, .(
    country,
    variant = get(variant_col)
  )])

  data.table(
    step = step,
    unique_variants = uniqueN(data[[variant_col]]),
    country_variant_pairs = nrow(country_variant_pairs),
    cases_sequences = sum(data$numerator, na.rm = TRUE)
  )
}

capture_command <- function(command, args) {
  output <- suppressWarnings(system2(command, args, stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    fail("Command failed: ", command, " ", paste(args, collapse = " "))
  }
  output
}

extract_version <- function(lines, label) {
  line <- lines[startsWith(tolower(lines), tolower(paste0(label, ":")))][1L]
  if (is.na(line)) return(NA_character_)
  trimws(sub("^[^:]+:", "", line))
}

first_token <- function(lineage) sub("\\..*$", "", lineage)

expand_one_alias <- function(lineage, alias_definitions, max_steps = 30L) {
  if (is.na(lineage) || !nzchar(lineage)) return(NA_character_)
  current <- lineage

  for (step in seq_len(max_steps)) {
    prefix <- first_token(current)
    if (prefix %in% c("A", "B")) return(current)
    definition <- alias_definitions[[prefix]]

    # Recombinant aliases have multiple parents and cannot be reversibly
    # expanded onto one tree path; retain their official X* hierarchy.
    if (is.null(definition) ||
        length(definition) != 1L ||
        !is.character(definition)) {
      return(current)
    }

    suffix <- substring(current, nchar(prefix) + 1L)
    next_value <- paste0(definition, suffix)
    if (identical(next_value, current)) return(current)
    current <- next_value
  }
  fail("Alias expansion exceeded ", max_steps, " steps for ", lineage)
}

lineage_depth <- function(lineage) {
  lengths(regmatches(lineage, gregexpr("\\.", lineage, fixed = FALSE)))
}

is_same_or_descendant <- function(lineage, root) {
  lineage == root | startsWith(lineage, paste0(root, "."))
}

is_terminal_in_country <- function(root, country_lineages) {
  !any(startsWith(country_lineages, paste0(root, ".")))
}

is_informative_root <- function(root) {
  is_recombinant_family <- grepl("^X[A-Z]", root)
  is_recombinant_family | lineage_depth(root) >= minimum_nonrecombinant_depth
}

coarsest_matching_root <- function(lineage, roots) {
  matches <- roots[is_same_or_descendant(lineage, roots)]
  if (length(matches) == 0L) return(lineage)
  matches[order(lineage_depth(matches), nchar(matches), matches)][1L]
}

match_who_direct <- function(lineage) {
  matches <- which(is_same_or_descendant(lineage, who_reference$canonical_root))
  if (length(matches) == 0L) return(NA_character_)
  best <- matches[which.max(nchar(who_reference$canonical_root[matches]))]
  who_reference$who_variant[best]
}

resolve_alias_parent_leaves <- function(alias_name,
                                        alias_definitions,
                                        seen = character(),
                                        max_steps = 30L) {
  if (max_steps <= 0L || alias_name %in% seen) return(alias_name)

  clean_name <- sub("\\*$", "", alias_name)
  prefix <- first_token(clean_name)
  definition <- alias_definitions[[prefix]]
  if (is.null(definition)) return(clean_name)

  if (is.character(definition) && length(definition) == 1L) {
    return(expand_one_alias(clean_name, alias_definitions))
  }
  if (is.list(definition) || length(definition) > 1L) {
    parents <- unlist(definition, use.names = FALSE)
    return(unlist(lapply(
      parents,
      resolve_alias_parent_leaves,
      alias_definitions = alias_definitions,
      seen = c(seen, alias_name),
      max_steps = max_steps - 1L
    )))
  }
  clean_name
}

match_who_recombinant <- function(lineage, alias_definitions) {
  if (!grepl("^X", lineage)) return(NA_character_)
  parent_leaves <- unique(resolve_alias_parent_leaves(
    first_token(lineage),
    alias_definitions
  ))
  parent_labels <- vapply(parent_leaves, match_who_direct, character(1))
  if (length(parent_leaves) > 0L &&
      all(!is.na(parent_labels)) &&
      uniqueN(parent_labels) == 1L) {
    return(unique(parent_labels))
  }
  NA_character_
}

# -----------------------------------------------------------------------------
# Validate the frozen Pango reference and save its exact snapshot
# -----------------------------------------------------------------------------

if (!file.exists(raw_data_file)) fail("Missing raw summary: ", raw_data_file)
if (!file.exists(country_file)) fail("Missing 15-country list: ", country_file)
if (!nzchar(pangolin_bin) || !file.exists(pangolin_bin)) {
  fail("Missing Pangolin executable. Add it to PATH or set PANGOLIN_BIN.")
}

all_versions <- capture_command(pangolin_bin, "--all-versions")
installed_pangolin <- extract_version(all_versions, "pangolin")
installed_data <- extract_version(all_versions, "pangolin-data")
if (!identical(installed_pangolin, expected_pangolin_version) ||
    !identical(installed_data, expected_pangolin_data_version)) {
  fail(
    "Frozen environment mismatch. Found Pangolin ", installed_pangolin,
    " and pangolin-data ", installed_data, "; expected ",
    expected_pangolin_version, " and ", expected_pangolin_data_version, "."
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(all_versions, file.path(output_dir, "pangolin_all_versions.txt"))

alias_text <- capture_command(pangolin_bin, "--aliases")
alias_snapshot_file <- file.path(output_dir, "alias_key_pangolin_data_1.39.json")
writeLines(alias_text, alias_snapshot_file)
alias_definitions <- read_json(alias_snapshot_file, simplifyVector = FALSE)

# -----------------------------------------------------------------------------
# Load all dates for the fixed 15-country set; do not use 01 start/end dates
# -----------------------------------------------------------------------------

country_rank <- fread(country_file, select = c("rank", "country"))
if (nrow(country_rank) != 15L || uniqueN(country_rank$country) != 15L) {
  fail("The country list must contain exactly 15 unique countries.")
}

raw_country_names <- copy(country_rank)
raw_country_names[country == "United States", raw_country := "USA"]
raw_country_names[is.na(raw_country), raw_country := country]

raw_data <- fread(
  raw_data_file,
  select = c("country", "date", "lineage", "denominator", "numerator")
)
analysis_data <- raw_data[country %in% raw_country_names$raw_country]
analysis_data[raw_country_names, on = .(country = raw_country), country_display := i.country]
analysis_data[, country := country_display]
analysis_data[, country_display := NULL]
analysis_data[, date := as.IDate(date)]

unassigned_summary <- analysis_data[
  !is.na(lineage) & tolower(trimws(lineage)) == "unassigned",
  .(
    excluded_rows = .N,
    excluded_sequences = sum(numerator, na.rm = TRUE)
  )
]
analysis_data <- analysis_data[
  is.na(lineage) | tolower(trimws(lineage)) != "unassigned"
]

missing_countries <- setdiff(country_rank$country, unique(analysis_data$country))
if (length(missing_countries) > 0L) {
  fail("No raw data found for: ", paste(missing_countries, collapse = ", "))
}

denominator_check <- unique(analysis_data[, .(country, date, denominator)])[
  , .N, by = .(country, date)
][N != 1L]
if (nrow(denominator_check) > 0L) {
  fail("Some country-date combinations have multiple denominators.")
}

# -----------------------------------------------------------------------------
# Construct one cross-country lineage mapping
# -----------------------------------------------------------------------------

lineage_mapping <- data.table(
  original_lineage = sort(unique(analysis_data$lineage))
)
lineage_mapping[, canonical_lineage := vapply(
  original_lineage,
  expand_one_alias,
  character(1),
  alias_definitions = alias_definitions
)]
lineage_mapping[, alias_expanded := canonical_lineage != original_lineage]

analysis_data[lineage_mapping, on = .(lineage = original_lineage),
              canonical_lineage := i.canonical_lineage]

country_lineages <- unique(analysis_data[, .(country, canonical_lineage)])
terminal_roots <- country_lineages[, {
  observed <- canonical_lineage
  .(
    canonical_root = observed[vapply(
      observed,
      is_terminal_in_country,
      logical(1),
      country_lineages = observed
    )]
  )
}, by = country]
terminal_roots[, informative := vapply(canonical_root, is_informative_root, logical(1))]

candidate_roots <- sort(unique(terminal_roots[informative == TRUE, canonical_root]))
has_candidate_ancestor <- vapply(candidate_roots, function(lineage) {
  other_roots <- setdiff(candidate_roots, lineage)
  any(is_same_or_descendant(lineage, other_roots))
}, logical(1))
global_common_roots <- candidate_roots[!has_candidate_ancestor]

lineage_mapping[, harmonized_group := vapply(
  canonical_lineage,
  coarsest_matching_root,
  character(1),
  roots = global_common_roots
)]
lineage_mapping[, was_merged := canonical_lineage != harmonized_group]
lineage_mapping[, levels_rolled_up := pmax(
  0L,
  lineage_depth(canonical_lineage) - lineage_depth(harmonized_group)
)]

analysis_data[lineage_mapping, on = .(lineage = original_lineage), `:=`(
  canonical_lineage = i.canonical_lineage,
  harmonized_group = i.harmonized_group
)]

# -----------------------------------------------------------------------------
# Aggregate all available dates and quantify changes for each country
# -----------------------------------------------------------------------------

daily_grouped <- analysis_data[, .(
  numerator = sum(numerator, na.rm = TRUE)
), by = .(country, date, denominator, harmonized_group)]
daily_grouped[, share := numerator / denominator]
setorder(daily_grouped, country, date, -numerator, harmonized_group)
if (any(daily_grouped$share > 1 + 1e-10, na.rm = TRUE)) {
  fail("A grouped daily share exceeds 1; aggregation is invalid.")
}

# Match harmonized groups to the 13 historical WHO variants. Direct
# ancestor/descendant matching is preferred; recombinant ancestry is used only
# when all resolved parents point to the same WHO group.
who_mapping <- data.table(
  harmonized_group = sort(unique(daily_grouped$harmonized_group))
)
who_mapping[, who_variant := vapply(
  harmonized_group,
  match_who_direct,
  character(1)
)]
who_mapping[is.na(who_variant) & grepl("^X", harmonized_group),
            who_variant := vapply(
              harmonized_group,
              match_who_recombinant,
              character(1),
              alias_definitions = alias_definitions
            )]

annotated_variants <- merge(
  daily_grouped,
  who_mapping,
  by = "harmonized_group",
  all.x = TRUE
)
setnames(annotated_variants, "harmonized_group", "variant")
setcolorder(
  annotated_variants,
  c("country", "date", "variant", "who_variant", "denominator", "numerator", "share")
)
setorder(annotated_variants, country, date, -numerator, variant)

if (nrow(annotated_variants) != nrow(daily_grouped) ||
    uniqueN(annotated_variants$country) != 15L) {
  fail("WHO annotation unexpectedly dropped rows or countries.")
}

# Apply one combined eligibility rule per country-variant trajectory:
#   peak share > 1% AND >=7 consecutive calendar days with share > 1%.
# Missing dates break a run because qualifying dates must be exactly one day
# apart. The two component flags are retained as diagnostics, but filtering is
# performed once using `meets_combined_rule`.
variant_periods <- annotated_variants[, .(
  peak_share = max(share, na.rm = TRUE),
  peak_share_date = date[which.max(share)][1L],
  first_observed_date = min(date),
  last_observed_date = max(date)
), by = .(country, variant, who_variant)]

qualifying_days <- annotated_variants[
  share > share_threshold,
  .(country, variant, who_variant, date)
]
setorder(qualifying_days, country, variant, date)
qualifying_days[, previous_date := shift(date), by = .(country, variant)]
qualifying_days[, run_id := cumsum(
  is.na(previous_date) | as.integer(date - previous_date) != 1L
), by = .(country, variant)]

qualifying_runs <- qualifying_days[, .(
  run_start = min(date),
  run_end = max(date),
  run_length_days = .N
), by = .(country, variant, who_variant, run_id)][
  run_length_days >= minimum_run_days
]

run_summary <- qualifying_runs[, .(
  first_7day_start = min(run_start),
  last_7day_end = max(run_end),
  n_qualifying_runs = .N
), by = .(country, variant, who_variant)]

variant_periods[run_summary, on = .(country, variant, who_variant), `:=`(
  first_7day_start = i.first_7day_start,
  last_7day_end = i.last_7day_end,
  n_qualifying_runs = i.n_qualifying_runs
)]
variant_periods[is.na(n_qualifying_runs), n_qualifying_runs := 0L]
variant_periods[, `:=`(
  peak_above_1pct = peak_share > share_threshold,
  has_confirmed_7day_period = n_qualifying_runs > 0L
)]
variant_periods[, meets_combined_rule :=
  peak_above_1pct & has_confirmed_7day_period]
variant_periods[, duration_days := fifelse(
  meets_combined_rule,
  as.integer(last_7day_end - first_7day_start) + 1L,
  NA_integer_
)]
variant_periods[, duration_over_30_days :=
  meets_combined_rule & duration_days > 30L]

# Keep only trajectories satisfying both parts of the single eligibility rule.
variant_periods <- variant_periods[meets_combined_rule == TRUE]
setorder(variant_periods, country, -peak_share, variant)

eligible_variants <- merge(
  annotated_variants,
  variant_periods,
  by = c("country", "variant", "who_variant"),
  all = FALSE
)
setcolorder(
  eligible_variants,
  c(
    "country", "date", "variant", "who_variant", "denominator", "numerator",
    "share", "peak_share", "peak_share_date", "first_7day_start",
    "last_7day_end", "duration_days", "has_confirmed_7day_period",
    "duration_over_30_days", "n_qualifying_runs", "first_observed_date",
    "last_observed_date", "peak_above_1pct", "meets_combined_rule"
  )
)
setorder(eligible_variants, country, variant, date)

if (uniqueN(eligible_variants[, .(country, variant)]) != nrow(variant_periods) ||
    any(!eligible_variants$meets_combined_rule)) {
  fail("Combined eligibility filtering failed validation.")
}

# Start every eligible trajectory on the first day of its earliest qualifying
# seven-day run. All later observations are retained, even if share later falls
# to 1% or below. Keep the same seven columns previously produced by step 03.
trimmed_variants <- eligible_variants[
  !is.na(first_7day_start) & date >= first_7day_start,
  .(
    country,
    date,
    variant,
    who_variant,
    denominator,
    numerator,
    share
  )
]
setorder(trimmed_variants, country, variant, date)

trimmed_starts <- trimmed_variants[, .(
  retained_start = min(date)
), by = .(country, variant)]
expected_starts <- unique(variant_periods[, .(
  country,
  variant,
  first_7day_start
)])
start_check <- trimmed_starts[expected_starts, on = .(country, variant)]

if (nrow(trimmed_starts) != nrow(variant_periods) ||
    any(start_check$retained_start != start_check$first_7day_start)) {
  fail("Trajectory start-date trimming failed validation.")
}

country_group_peaks <- daily_grouped[, .(
  total_sequences = sum(numerator),
  first_observed_date = min(date),
  last_observed_date = max(date),
  peak_share = max(share),
  peak_share_date = date[which.max(share)][1L]
), by = .(country, harmonized_group)]
country_group_peaks[, `:=`(
  ever_above_1pct = peak_share > share_threshold,
  ever_above_20pct = peak_share > diagnostic_high_share_threshold
)]

country_counts <- analysis_data[, .(
  first_date = min(date),
  last_date = max(date),
  n_raw_lineages = uniqueN(lineage),
  n_alias_normalized_lineages = uniqueN(canonical_lineage),
  n_harmonized_groups = uniqueN(harmonized_group),
  n_raw_labels_merged = uniqueN(lineage) - uniqueN(harmonized_group)
), by = country]
country_counts[country_group_peaks[, .(
  n_groups_ever_above_1pct = sum(ever_above_1pct),
  n_groups_ever_above_20pct = sum(ever_above_20pct)
), by = country], on = "country", `:=`(
  n_groups_ever_above_1pct = i.n_groups_ever_above_1pct,
  n_groups_ever_above_20pct = i.n_groups_ever_above_20pct
)]
country_counts[, reduction_percent :=
  100 * n_raw_labels_merged / n_raw_lineages]
country_counts[country_rank, on = "country", rank := i.rank]
setorder(country_counts, rank)

eligible_country_counts <- variant_periods[, .(
  n_eligible_variants = .N,
  n_variants_duration_over_30_days = sum(duration_over_30_days),
  median_duration_days = as.numeric(median(duration_days, na.rm = TRUE)),
  mean_duration_days = mean(duration_days, na.rm = TRUE),
  maximum_duration_days = max(duration_days, na.rm = TRUE),
  n_who_variants_present = uniqueN(who_variant, na.rm = TRUE),
  n_unmatched_eligible_variants = sum(is.na(who_variant))
), by = country]
eligible_country_counts[country_counts, on = "country", `:=`(
  n_harmonized_groups_before_eligibility_filter = i.n_harmonized_groups,
  percent_of_harmonized_groups_eligible =
    100 * n_eligible_variants / i.n_harmonized_groups,
  rank = i.rank
)]
setorder(eligible_country_counts, rank)

root_summary <- terminal_roots[canonical_root %in% global_common_roots,
  .(
    n_countries_where_terminal = uniqueN(country),
    countries_where_terminal = paste(sort(unique(country)), collapse = "; ")
  ),
  by = canonical_root
][order(lineage_depth(canonical_root), canonical_root)]

mapping_country_presence <- unique(analysis_data[, .(
  country,
  original_lineage = lineage,
  canonical_lineage,
  harmonized_group
)])
mapping_country_presence[lineage_mapping, on = .(original_lineage), `:=`(
  alias_expanded = i.alias_expanded,
  was_merged = i.was_merged,
  levels_rolled_up = i.levels_rolled_up
)]
setorder(mapping_country_presence, country, harmonized_group, canonical_lineage)

# -----------------------------------------------------------------------------
# Save data and a compact run manifest
# -----------------------------------------------------------------------------

fwrite(country_counts, file.path(output_dir, "country_counts_before_after.csv"))
fwrite(lineage_mapping, file.path(output_dir, "global_lineage_mapping.csv"))
fwrite(
  mapping_country_presence,
  file.path(output_dir, "country_lineage_mapping.csv")
)
fwrite(terminal_roots, file.path(output_dir, "country_terminal_roots.csv"))
fwrite(root_summary, file.path(output_dir, "global_common_roots.csv"))
fwrite(daily_grouped, file.path(output_dir, "country_date_grouped_counts.csv"))
fwrite(country_group_peaks, file.path(output_dir, "country_group_peak_share.csv"))
fwrite(who_reference, file.path(output_dir, "who_roots.csv"))
fwrite(who_mapping, file.path(output_dir, "who_mapping.csv"))
fwrite(
  eligible_country_counts,
  file.path(output_dir, "country_eligibility_counts.csv")
)
fwrite(variant_periods, file.path(output_dir, "variant_periods.csv"))
saveRDS(trimmed_variants, output_rds_file, compress = TRUE)

manifest <- list(
  created_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  method = paste(
    "Alias expansion using frozen pangolin-data; country-terminal reporting",
    "roots; coarsest terminal root shared globally along each lineage chain"
  ),
  warning = paste(
    "This is label harmonization only. It cannot determine whether an ancestral",
    "label represents a true ancestral sequence or an unresolved descendant."
  ),
  uses_01_start_end = FALSE,
  raw_data_file = raw_data_file,
  raw_data_md5 = unname(tools::md5sum(raw_data_file)),
  country_file_used_for_names_only = country_file,
  pangolin_version = installed_pangolin,
  pangolin_data_version = installed_data,
  alias_snapshot_md5 = unname(tools::md5sum(alias_snapshot_file)),
  minimum_nonrecombinant_depth = minimum_nonrecombinant_depth,
  who_annotation = paste(
    "13 historical WHO variants represented by 14 canonical roots; descendants",
    "matched; X* matched only when all resolved parents have one WHO label;",
    "unmatched variants retained with who_variant = NA"
  ),
  primary_rds_file = output_rds_file,
  eligibility_rule = paste(
    "country-specific peak daily share strictly greater than 0.01 AND at least",
    "one run of seven consecutive calendar days with share strictly greater",
    "than 0.01"
  ),
  trajectory_start_rule = paste(
    "retain observations from the first day of the earliest qualifying",
    "seven-day run onward"
  ),
  period_definition = paste(
    "First day of earliest >=7-calendar-day run with share >0.01 through",
    "last day of latest qualifying run"
  ),
  long_duration_definition = "duration_days > 90",
  countries = country_counts$country,
  n_global_common_roots = length(global_common_roots)
)
write_json(
  manifest,
  file.path(output_dir, "run_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

# Print only the compact, consistently defined retention summary to the console.
step_summary <- rbindlist(list(
  summarize_step("1. Selected 15-country data", analysis_data, "lineage"),
  summarize_step("2. Harmonized variants", daily_grouped, "harmonized_group"),
  summarize_step("3. WHO annotation (no filtering)", annotated_variants, "variant"),
  summarize_step(
    "4. Combined eligibility rule",
    eligible_variants,
    "variant"
  ),
  summarize_step(
    "5. Trim before first qualifying run",
    trimmed_variants,
    "variant"
  )
))

# Compare every step with the immediately preceding step. Negative changes mean
# that labels, country-variant pairs, or sequences decreased at that step.
step_summary[, `:=`(
  unique_variants_change = unique_variants - shift(unique_variants),
  country_variant_pairs_change =
    country_variant_pairs - shift(country_variant_pairs),
  cases_sequences_change = cases_sequences - shift(cases_sequences),
  cases_sequences_retained_pct = round(
    100 * cases_sequences / shift(cases_sequences),
    2
  )
)]
setcolorder(step_summary, c(
  "step",
  "unique_variants", "unique_variants_change",
  "country_variant_pairs", "country_variant_pairs_change",
  "cases_sequences", "cases_sequences_change",
  "cases_sequences_retained_pct"
))

cat("\nProcessing summary (change = current step - previous step)\n")
print(step_summary)

message("Saved trimmed eligible variant trajectories: ", output_rds_file)
message(
  "Excluded Unassigned records before harmonization: ",
  format(unassigned_summary$excluded_rows, big.mark = ","), " rows; ",
  format(unassigned_summary$excluded_sequences, big.mark = ","), " sequences"
)
message("Saved diagnostics under: ", output_dir)
