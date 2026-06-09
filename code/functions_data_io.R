#' Prepare early-period wide matrix of lineage share by country
#'
#' Builds a country × day feature table from long time-series, aligned to the
#' first day the lineage exceeds a minimum share threshold. Keeps countries with
#' at least `input_days` of data and where the lineage never becomes dominant
#' (max share < `dominant_threshold`) within that window.
#'
#' @param data A data.frame/tibble with columns: country_label, date, share.
#' @param input_days Integer >= 1. Number of days to keep starting from first detection.
#' @param dominant_threshold Numeric in (0, 1]. Max share must be < this to keep a row.
#' @param min_share Numeric in (0, 1]. Threshold to define first detection (default 0.01).
#' @param label_sep Separator used to split country_label into country and classified_label.
#' @param values_fill Value used to fill missing days in pivot (default NA_real_).
#'
#' @return A tibble with columns: country, classified_label, day_1 ... day_{input_days}.
#' @export
# Prepare early-period wide matrix of lineage share by country
prepare_model_data <- function(data,
                               input_days,
                               dominant_threshold = 0.5,
                               min_share = 0.01,
                               values_fill = NA_real_) {
  req <- c("country","lineage","date","share")
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  if (!is.numeric(input_days) || length(input_days) != 1L || input_days < 1)
    stop("input_days must be a single integer >= 1")
  if (!inherits(data$date, "Date")) {
    data$date <- as.Date(data$date)
    if (any(is.na(data$date))) stop("`date` must be coercible to Date")
  }
  
  data_clean <- data %>%
    dplyr::group_by(country, lineage, date) %>%
    dplyr::summarise(share = mean(share, na.rm = TRUE), .groups = "drop")
  
  day_names <- paste0("day_", seq_len(input_days))
  combos <- dplyr::distinct(data_clean, country, lineage)
  
  first_tbl <- data_clean %>%
    dplyr::group_by(country, lineage) %>%
    dplyr::summarise(
      any_above = any(share >= min_share, na.rm = TRUE),
      first_share_date = if (any(share >= min_share, na.rm = TRUE))
        min(date[share >= min_share], na.rm = TRUE) else as.Date(NA),
      .groups = "drop"
    )
  
  to_exclude <- first_tbl %>% dplyr::filter(!any_above) %>%
    dplyr::transmute(key = paste(country, lineage, sep = " | ")) %>% dplyr::pull()
  if (length(to_exclude) > 0)
    warning(sprintf("No observation reached %.3f for %d combos (e.g., %s ...)",
                    min_share, length(to_exclude), paste(utils::head(to_exclude, 5L), collapse = "; ")))
  
  data_kept <- dplyr::inner_join(
    data_clean,
    dplyr::filter(first_tbl, any_above),
    by = c("country","lineage")
  )
  
  out_good <- data_kept %>%
    dplyr::mutate(time = as.integer(date - first_share_date) + 1L) %>%
    dplyr::filter(!is.na(time), time > 0L, time <= input_days) %>%
    dplyr::group_by(country, lineage) %>%
    dplyr::filter(max(time, na.rm = TRUE) >= input_days,
                  max(share, na.rm = TRUE) < dominant_threshold) %>%
    dplyr::ungroup() %>%
    dplyr::select(country, lineage, time, share) %>%
    tidyr::pivot_wider(names_from = time, values_from = share,
                       names_prefix = "day_", values_fill = values_fill) %>%
    dplyr::distinct(country, lineage, .keep_all = TRUE)
  
  missing_days <- setdiff(day_names, names(out_good))
  for (nm in missing_days) out_good[[nm]] <- NA_real_
  
  X <- dplyr::left_join(combos, out_good, by = c("country","lineage")) %>%
    dplyr::relocate(dplyr::all_of(c("country","lineage", day_names))) %>%
    tibble::as_tibble()
  
  window_eval <- data_kept %>%
    dplyr::mutate(time = as.integer(date - first_share_date) + 1L) %>%
    dplyr::filter(!is.na(time), time > 0L, time <= input_days) %>%
    dplyr::group_by(country, lineage) %>%
    dplyr::summarise(
      window_len_ge_input = max(time, na.rm = TRUE) >= input_days,
      max_share_window    = max(share, na.rm = TRUE),
      reached_dominant_in_window = max_share_window >= dominant_threshold,
      .groups = "drop"
    )
  
  combo_flags <- combos %>%
    dplyr::left_join(first_tbl, by = c("country","lineage")) %>%
    dplyr::left_join(window_eval, by = c("country","lineage")) %>%
    dplyr::mutate(
      window_len_ge_input = dplyr::coalesce(window_len_ge_input, FALSE),
      reached_dominant_in_window = dplyr::coalesce(reached_dominant_in_window, FALSE),
      max_share_window = ifelse(is.finite(max_share_window), max_share_window, NA_real_)
    ) %>%
    dplyr::rename(any_above_min_share = any_above)
  
  list(
    X = X,
    combo_flags = combo_flags
  )
}





#' Extract growth/trend features from early-period variant share table
#'
#' Computes a rich set of handcrafted growth- and trend-related features from a
#' wide country × day table (e.g., \code{model_data_14}), including slopes,
#' curvature, volatility, peaks, area-under-curve, time-to-thresholds, ACF/PACF,
#' monotonic run length, sign changes, and doubling time. Optionally augments
#' with automated features from \pkg{Rcatch22}/\pkg{catch22} (the 22 canonical
#' time-series features).
#'
#' @param model_data A data.frame/tibble with columns \code{country},
#'   \code{classified_label}, and day-wise numeric columns (e.g., \code{day_1}, \code{day_2}, ...).
#' @param day_prefix Regular expression to identify day columns (default \code{"^day_\\\\d+$"}).
#' @param eps Small numeric constant to avoid \code{log(0)} or division-by-zero in transforms.
#' @param early_ks Integer vector of window lengths for early-window slope calculations.
#' @param thresholds Numeric vector of share thresholds (in proportion, e.g. 0.05 for 5\%)
#'   for which to compute "time-to-threshold" features.
#' @param add_auto_features Logical; if \code{TRUE}, compute 22 features using
#'   \pkg{Rcatch22} (preferred if installed) or fallback to \pkg{catch22}.
#'
#' @return A tibble keyed by \code{country} and \code{classified_label} containing
#'   handcrafted feature columns (numeric) and, if requested, any available
#'   \pkg{catch22} features (prefixed with \code{c22_}).
#' @export
extract_features_from_model_data <- function(model_data,
                                             day_prefix        = "^day_\\d+$",
                                             eps               = 1e-6,
                                             early_ks          = c(7),
                                             add_auto_features = TRUE) {
  stopifnot(all(c("country","lineage") %in% names(model_data)))
  day_cols <- grep(day_prefix, names(model_data), value = TRUE)
  if (length(day_cols) < 5) stop("not enough day_* columns found")
  day_cols <- day_cols[order(as.integer(sub("^day_", "", day_cols)))]
  
  .lin_slope <- function(x, t) {
    if (all(!is.finite(x))) return(NA_real_)
    as.numeric(tryCatch(stats::coef(stats::lm(x ~ t))[2], error = function(e) NA_real_))
  }
  .win_slope <- function(x, k) {
    k <- min(k, length(x)); if (k < 2) return(NA_real_)
    t <- seq_len(k); .lin_slope(x[t], t)
  }
  .longest_inc_run <- function(x) {
    d <- diff(x); s <- as.integer(d > 0L)
    if (!any(s == 1L, na.rm = TRUE)) return(0L)
    r <- rle(s); max(c(0L, r$lengths[r$values == 1L]))
  }
  
  keys <- model_data[, c("country","lineage")]
  X    <- as.data.frame(model_data[, day_cols, drop = FALSE])
  
  all_na_rows <- apply(X, 1, function(r) all(is.na(r)))
  n_skip <- sum(all_na_rows)
  if (n_skip > 0) warning(sprintf("Skipped %d rows with all-NA day_* values.", n_skip))
  if (n_skip == nrow(X)) return(tibble::tibble())
  
  keys <- keys[!all_na_rows, , drop = FALSE]
  X    <- X[!all_na_rows, , drop = FALSE]
  
  feat_list <- lapply(seq_len(nrow(X)), function(i) {
    v0 <- as.numeric(X[i, ])
    zero_count <- sum(v0 == 0, na.rm = TRUE)
    
    v <- v0
    v[!is.finite(v) | is.na(v)] <- 0
    t <- seq_along(v)
    
    v_log <- log(pmax(v, 0) + eps)
    
    slope_log <- .lin_slope(v_log, t)
    win_slopes <- setNames(
      sapply(early_ks, function(k) .win_slope(v, k)),
      paste0("slope_y_1_", pmin(early_ks, length(v)))
    )
    
    peak_val <- max(v, na.rm = TRUE)
    auc_norm <- sum(v, na.rm = TRUE) / length(v)
    longest_inc_run <- .longest_inc_run(v)
    
    c(
      slope_log         = slope_log,
      win_slopes,
      peak_val          = peak_val,
      auc_norm          = auc_norm,
      longest_inc_run   = longest_inc_run,
      zero_count        = zero_count
    )
  })
  
  handcrafted <- as.data.frame(do.call(rbind, feat_list), stringsAsFactors = FALSE)
  handcrafted[] <- lapply(handcrafted, function(col) if (is.character(col)) suppressWarnings(as.numeric(col)) else col)
  out <- cbind(keys, handcrafted, row.names = NULL)
  
  if (isTRUE(add_auto_features)) {
    has_Rcatch22 <- requireNamespace("Rcatch22", quietly = TRUE)
    has_catch22  <- requireNamespace("catch22",  quietly = TRUE)
    if (has_Rcatch22 || has_catch22) {
      .c22 <- function(x) {
        x[!is.finite(x) | is.na(x)] <- 0
        fn <- if (has_Rcatch22) get("catch22_all", envir = asNamespace("Rcatch22"))
        else               get("catch22_all", envir = asNamespace("catch22"))
        fmls <- names(formals(fn)); args <- list()
        if (length(fmls) >= 1L) args[[ fmls[1] ]] <- x else args[[1]] <- x
        if ("catch24"   %in% fmls) args$catch24   <- FALSE
        if ("normalize" %in% fmls) args$normalize <- FALSE
        if ("silent"    %in% fmls) args$silent    <- TRUE
        if ("doParallel"%in% fmls) args$doParallel<- FALSE
        res <- tryCatch(do.call(fn, args), error = function(e) NULL)
        if (is.null(res)) return(setNames(as.numeric(NA), character(0)))
        if (is.data.frame(res)) {
          val_col  <- grep("value", names(res),  ignore.case = TRUE, value = TRUE)[1]
          name_col <- grep("name",  names(res),  ignore.case = TRUE, value = TRUE)[1]
          if (!is.na(val_col) && !is.na(name_col)) {
            v <- res[[val_col]]; nm <- res[[name_col]]; names(v) <- nm
          } else { v <- res[[2]]; nm <- res[[1]]; names(v) <- nm }
        } else if (is.numeric(res) && !is.null(names(res))) {
          v <- res
        } else {
          v <- as.numeric(res); names(v) <- paste0("c22_", seq_along(v))
        }
        names(v) <- paste0("c22_", make.names(names(v), unique = TRUE))
        v
      }
      c22_list  <- lapply(seq_len(nrow(X)), function(i) .c22(as.numeric(X[i, ])))
      all_names <- unique(unlist(lapply(c22_list, names)))
      c22_mat   <- do.call(rbind, lapply(c22_list, function(v) {
        outv <- rep(NA_real_, length(all_names)); names(outv) <- all_names
        outv[names(v)] <- as.numeric(v); outv
      }))
      out <- cbind(out, as.data.frame(c22_mat), row.names = NULL)
    } else {
      warning("Neither 'Rcatch22' nor 'catch22' is installed; skipping catch22 features.")
    }
  }
  
  tibble::as_tibble(out)
}



#' Clean `country_label` into `country` and `classified_label`
#'
#' Splits the \code{country_label} column (formatted like "Country-Variant")
#' into separate \code{country} and \code{classified_label} columns. If the
#' \code{classified_label} contains a Pango lineage in parentheses of the form
#' \code{(B.x.y...)} it is kept as-is; otherwise, any parentheses and contents
#' are removed.
#'
#' @param df A data.frame/tibble containing a \code{country_label} column.
#' @param sep Character. Separator between country and label (default "-").
#'
#' @return A tibble with new columns \code{country} and cleaned \code{classified_label}.
#' @export
clean_classified_label <- function(df, sep = "-") {
  if (!"country_label" %in% names(df)) {
    stop("`df` must contain a column named 'country_label'.")
  }
  
  df %>%
    tidyr::separate(country_label,
                    into = c("country", "classified_label"),
                    sep = sep) %>%
    dplyr::mutate(
      classified_label = ifelse(
        stringr::str_detect(.data$classified_label,
                            "\\(B\\.[0-9]+(\\.[0-9]+)*\\)"),
        .data$classified_label,
        stringr::str_remove(.data$classified_label, "\\s*\\(.*\\)")
      )
    )
}


#' 
#' #' Classify lineages to variant labels and aggregate daily shares (per country)
#' #'
#' #' Given a country's lineage-level time series, this function maps PANGO
#' #' lineages to higher-level variant labels (using an alias list and a helper
#' #' \code{get_classification()}), filters to sufficiently supported lineages,
#' #' and returns a date × classified\_label table with aggregated counts and
#' #' shares. It also returns the classification table used for the mapping.
#' #'
#' #' The mapping is performed via a prefix (starts-with) fuzzy join on
#' #' \emph{decimalized} lineage strings (a trailing dot is appended) so that,
#' #' for example, \code{"BA.2.86.1"} correctly inherits the classification for
#' #' \code{"BA.2.86"}. When multiple candidate matches exist, the longest
#' #' matching prefix is retained. To avoid over-aggregation, we only accept
#' #' mappings where the lineage depth differs from the classified label by at
#' #' most one level; otherwise the lineage is set to \code{Other} and dropped.
#' #'
#' #' @param data_cty A data frame (single country) with columns:
#' #'   \describe{
#' #'     \item{\code{country}}{Country name (character).}
#' #'     \item{\code{date}}{Sampling date (Date or coercible).}
#' #'     \item{\code{lineage}}{PANGO lineage (character), e.g., \code{"BA.2.86.1"}.}
#' #'     \item{\code{numerator}}{Number of sequences for this lineage/date (integer).}
#' #'     \item{\code{denominator}}{Total sequences for this date (integer).}
#' #'     \item{\code{share}}{Lineage share for this date (\eqn{\in[0,1]}).}
#' #'   }
#' #' @param alias_list A lookup object consumed by \code{get_classification()}
#' #'   that defines how raw lineages are grouped to \code{classified_label} and
#' #'   \code{classified_unasliased}. Its exact structure should match the helper.
#' #' @param p_lim Numeric scalar in \eqn{[0,1]}. Minimum global proportion
#' #'   (across all dates) for a lineage to be retained in the classification
#' #'   table produced by \code{get_classification()}. Default \code{0.001}.
#' #' @param n_lim Integer scalar. Minimum total sequence count for a lineage to be
#' #'   retained in the classification table produced by \code{get_classification()}.
#' #'   Default \code{50}.
#' #'
#' #' @return A named list with two elements:
#' #' \describe{
#' #'   \item{\code{data}}{A tibble with columns
#' #'     \code{date}, \code{classified_label}, \code{numerator}, \code{share},
#' #'     \code{denominator}; aggregated per date × classified\_label.}
#' #'   \item{\code{classification}}{The classification table returned by
#' #'     \code{get_classification()} (augmented with \code{decimal_lineage}).}
#' #' }
#' #'
#' #' @details
#' #' The function:
#' #' \enumerate{
#' #'   \item Summarizes total counts per lineage to inform \code{get_classification()}.
#' #'   \item Builds a classification table under thresholds \code{p_lim}, \code{n_lim}.
#' #'   \item Performs a prefix fuzzy join on decimalized lineages to assign
#' #'   \code{classified_label}.
#' #'   \item Keeps the longest matching classification per raw lineage.
#' #'   \item Ensures lineage depth differs by at most one level from the assigned
#' #'   classification; otherwise labeled \code{Other} and removed.
#' #'   \item Aggregates daily \code{numerator} and \code{share} over the same
#' #'   \code{date} × \code{classified_label}.
#' #' }
#' #'
#' #' Note: \code{share} is summed as provided. If your \code{share} is defined as
#' #' \code{numerator/denominator} at lineage level, summing shares across lineages
#' #' that share a \code{classified_label} is equivalent to recomputing
#' #' \code{sum(numerator) / denominator} when the denominator is common per date.
#' #' If denominators can vary within a date group, recompute \code{share}
#' #' after aggregation.
#' #'
#' #' @examples
#' #' \dontrun{
#' #' out <- classify_and_summarise_variants(
#' #'   data_cty = df_country,
#' #'   alias_list = alias_map,
#' #'   p_lim = 0.001,
#' #'   n_lim = 50
#' #' )
#' #' head(out$data)
#' #' head(out$classification)
#' #' }
#' #'
#' #' @importFrom dplyr group_by summarise mutate filter distinct select arrange
#' #' @importFrom dplyr if_else first
#' #' @importFrom stringr str_starts str_length str_count
#' #' @importFrom fuzzyjoin fuzzy_left_join
#' #' @export
#' classify_and_summarise_variants <- function(data_cty, alias_list, p_lim = 0.001, n_lim = 50) {
#'   summarised <- data_cty %>%
#'     dplyr::group_by(lineage) %>%
#'     dplyr::summarise(n = sum(numerator), .groups = "drop")
#'   
#'   classification_to_append <- get_classification(
#'     lineages = summarised$lineage,
#'     number_sequences = summarised$n,
#'     p_lim = p_lim,
#'     n_lim = n_lim,
#'     alias_list = alias_list
#'   ) %>%
#'     dplyr::mutate(decimal_lineage = paste0(classified_unasliased, "."))
#'   
#'   data_cty_classfied <- data_cty %>%
#'     dplyr::mutate(decimal_lineage = paste0(lineage, ".")) %>%
#'     fuzzyjoin::fuzzy_left_join(
#'       classification_to_append,
#'       by = "decimal_lineage",
#'       match_fun = stringr::str_starts
#'     ) %>%
#'     dplyr::filter(
#'       is.na(classified_unasliased) |
#'         stringr::str_length(classified_unasliased) ==
#'         max(stringr::str_length(classified_unasliased)),
#'       .by = lineage
#'     ) %>%
#'     dplyr::distinct(lineage, date, .keep_all = TRUE) %>%
#'     dplyr::mutate(
#'       length_class = stringr::str_count(classified_unasliased, "\\."),
#'       length_lineage = stringr::str_count(lineage, "\\."),
#'       classified_label = dplyr::if_else(
#'         (length_lineage - length_class) <= 1,
#'         classified_label,
#'         NA_character_
#'       ),
#'       classified_label = dplyr::if_else(
#'         is.na(classified_label), "Other", classified_label
#'       )
#'     ) %>%
#'     dplyr::filter(classified_label != "Other") %>%
#'     dplyr::select(country, date, lineage, numerator, denominator, share, classified_label)
#'   
#'   data_cty_classfied <- data_cty_classfied %>%
#'     dplyr::group_by(date, classified_label) %>%
#'     dplyr::summarise(
#'       numerator = sum(numerator, na.rm = TRUE),
#'       share = sum(share, na.rm = TRUE),   # recompute if your denominator varies within-day
#'       denominator = dplyr::first(denominator),
#'       .groups = "drop"
#'     )
#'   
#'   list(
#'     data = data_cty_classfied,
#'     classification = classification_to_append
#'   )
#' }
#' 
#' 
#' #' Classify, expand, and flatten variant time series across countries
#' #'
#' #' This function processes a full dataset of variant sequencing data by:
#' #' \enumerate{
#' #'   \item Iterating over each country
#' #'   \item Calling \code{classify_and_summarise_variants()} for lineage grouping
#' #'   \item Expanding time series so that each variant label has continuous
#' #'         daily coverage (missing dates filled with zeros)
#' #'   \item Flattening all country-variant series into a single tidy table
#' #' }
#' #'
#' #' @param data A data frame containing sequencing records with columns:
#' #'   \code{country}, \code{date}, \code{lineage}, \code{numerator},
#' #'   \code{denominator}, and \code{share}.
#' #' @param alias_list Lineage alias mapping passed to
#' #'   \code{classify_and_summarise_variants()}.
#' #' @param p_lim Numeric. Minimum global proportion for a lineage
#' #'   to be considered (default \code{0.001}).
#' #' @param n_lim Integer. Minimum number of sequences for a lineage
#' #'   to be considered (default \code{50}).
#' #'
#' #' @return A tibble with columns:
#' #'   \itemize{
#' #'     \item{\code{country_label}}{Concatenated country and classified label}
#' #'     \item{\code{date}}{Daily date sequence}
#' #'     \item{\code{numerator}}{Number of variant sequences (0 if missing)}
#' #'     \item{\code{denominator}}{Total sequences that day}
#' #'     \item{\code{share}}{Proportion of sequences belonging to the variant}
#' #'   }
#' #'
#' #' @examples
#' #' \dontrun{
#' #' classified_data <- process_classified_variants(
#' #'   data = data3,
#' #'   alias_list = alias_map
#' #' )
#' #' head(classified_data)
#' #' }
#' #'
#' #' @export
#' process_classified_variants <- function(data, alias_list,
#'                                         p_lim = 0.001, n_lim = 50) {
#'   country_list <- unique(data$country)
#'   classified_result_list <- list()
#'   
#'   # Step 1: classify per country
#'   for (cty in country_list) {
#'     message("Processing: ", cty)
#'     data_cty <- data %>% dplyr::filter(country == cty)
#'     classified_result_list[[cty]] <- classify_and_summarise_variants(
#'       data_cty, alias_list, p_lim = p_lim, n_lim = n_lim
#'     )
#'   }
#'   
#'   # Step 2: expand daily time series per country/label
#'   nested_classified_list <- list()
#'   for (cty in country_list) {
#'     df <- classified_result_list[[cty]]$data
#'     denom_by_date <- df %>% dplyr::distinct(date, denominator)
#'     
#'     nested_classified_list[[cty]] <- list()
#'     label_list <- unique(df$classified_label)
#'     
#'     for (label in label_list) {
#'       df_label <- df %>%
#'         dplyr::filter(classified_label == label) %>%
#'         dplyr::arrange(date)
#'       
#'       if (nrow(df_label) == 0) next
#'       
#'       date_range <- range(df_label$date, na.rm = TRUE)
#'       all_dates <- seq.Date(from = date_range[1], to = date_range[2], by = "day")
#'       
#'       complete_df <- tidyr::expand_grid(date = all_dates, classified_label = label) %>%
#'         dplyr::left_join(denom_by_date, by = "date") %>%
#'         dplyr::left_join(df_label, by = c("date", "classified_label")) %>%
#'         dplyr::mutate(
#'           country = cty,
#'           numerator = ifelse(is.na(numerator), 0, numerator),
#'           share = ifelse(is.na(share), 0, share),
#'           denominator = dplyr::coalesce(denominator.y, denominator.x)
#'         ) %>%
#'         dplyr::select(country, date, classified_label, numerator, denominator, share) %>%
#'         dplyr::arrange(date)
#'       
#'       nested_classified_list[[cty]][[label]] <- complete_df
#'     }
#'   }
#'   
#'   # Step 3: flatten nested list into one big tibble
#'   flat_df_list <- list()
#'   for (cty in names(nested_classified_list)) {
#'     for (label in names(nested_classified_list[[cty]])) {
#'       df <- nested_classified_list[[cty]][[label]]
#'       flat_df_list[[paste(cty, label, sep = "|")]] <- df
#'     }
#'   }
#'   
#'   classified_data <- dplyr::bind_rows(flat_df_list) %>%
#'     dplyr::mutate(country_label = paste0(country, "-", classified_label)) %>%
#'     dplyr::select(country_label, date, numerator, denominator, share)
#'   
#'   return(classified_data)
#' }



#' Plot smoothed stacked variant shares with a top-N legend
#'
#' @description
#' Given a single-country lineage time series (`sample_data`), this function
#' (i) completes missing date × label combinations, (ii) applies a centered
#' k-day moving average to shares, (iii) re-normalizes shares per day to sum
#' to 1, and (iv) draws a stacked area chart for all variants while showing
#' only the top-N variants (by weekly cumulative share) in the legend.
#'
#' @param sample_data A data frame for **one country** with columns:
#'   \code{country}, \code{classified_label}, \code{date},
#'   \code{numerator}, \code{denominator}, \code{share}.
#'   \code{date} must be coercible to \code{Date}.
#' @param cols A named character vector of colors where
#'   \code{names(cols)} are variant labels; any missing labels will be filled
#'   with distinct colors automatically.
#' @param top_n Integer. Number of variants to keep in the legend (default \code{10}).
#' @param k Integer. Window size (days) for centered moving average (default \code{7}).
#' @param legend_ncol Integer. Number of columns in the legend (default \code{1}).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' # Assume `sample_data` is filtered to a single country and contains required cols
#' # and `cols` is a named color vector for variants:
#' p <- plot_variant_stack_smoothed(sample_data, cols, top_n = 20, k = 7, legend_ncol = 1)
#' print(p)
#' }
#'
#' @importFrom dplyr mutate filter select group_by summarise arrange pull if_else
#' @importFrom dplyr ungroup left_join desc slice_head
#' @importFrom tidyr complete
#' @importFrom zoo rollmean
#' @importFrom ggplot2 ggplot aes geom_area scale_y_continuous scale_x_date
#' @importFrom ggplot2 scale_fill_manual labs theme_minimal theme element_blank
#' @importFrom ggplot2 guide_legend
#' @importFrom scales percent_format
#' @export
plot_variant_stack_smoothed <- function(sample_data,
                                        cols,
                                        legend_ncol = 1,
                                        k = 7,
                                        start_date,
                                        end_date,
                                        legends,
                                        legend_title = NULL) {
  stopifnot(all(c("country","lineage","date","numerator","denominator","share") %in% names(sample_data)))
  sample_data <- sample_data %>% dplyr::mutate(date = as.Date(date))
  cty <- unique(sample_data$country)
  if (length(cty) != 1L) stop("`sample_data` must contain exactly one country.")
  
  sample_data <- sample_data %>%
    dplyr::group_by(date) %>%
    dplyr::filter(sum(denominator, na.rm = TRUE) > 0 | sum(share, na.rm = TRUE) > 0) %>%
    dplyr::ungroup()
  
  all_dates   <- seq(min(sample_data$date), max(sample_data$date), by = "day")
  lineage_all <- sort(unique(sample_data$lineage))
  
  base_df <- sample_data %>%
    tidyr::complete(date = all_dates, lineage = lineage_all,
                    fill = list(numerator = 0, denominator = 0, share = 0)) %>%
    dplyr::arrange(lineage, date)
  
  smooth_df <- base_df %>%
    dplyr::group_by(lineage) %>%
    dplyr::mutate(
      share_smooth = zoo::rollapply(share, width = k, FUN = mean,
                                    align = "center", partial = TRUE, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(share_smooth = dplyr::coalesce(share_smooth, 0),
                  share_smooth = pmax(share_smooth, 0))
  
  nonothers_adj <- smooth_df %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(sum_s = sum(share_smooth, na.rm = TRUE),
                  fac   = ifelse(sum_s > 1, 1/sum_s, 1),
                  share = share_smooth * fac) %>%
    dplyr::ungroup() %>%
    dplyr::select(date, lineage, share)
  
  others_df <- nonothers_adj %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(sum_adj = sum(share, na.rm = TRUE), .groups = "drop") %>%
    dplyr::transmute(date, lineage = "others", share = pmax(0, 1 - sum_adj))
  
  plot_df <- dplyr::bind_rows(nonothers_adj, others_df)
  
  ord <- plot_df %>%
    dplyr::group_by(lineage) %>%
    dplyr::summarise(t_peak = date[which.max(share)], .groups = "drop") %>%
    dplyr::arrange(t_peak) %>%
    dplyr::pull(lineage)
  ord <- c(setdiff(ord, "others"), "others")
  plot_df <- plot_df %>% dplyr::mutate(lineage = factor(lineage, levels = ord))
  
  if (!"others" %in% names(cols)) cols <- c(cols, others = "lightgray") else cols["others"] <- "lightgray"
  miss_cols <- setdiff(levels(plot_df$lineage), names(cols))
  if (length(miss_cols)) cols <- c(cols, setNames(rep("grey80", length(miss_cols)), miss_cols))
  
  legends <- as.character(legends)
  present_lineages <- levels(plot_df$lineage)
  not_in_data <- setdiff(legends, present_lineages)
  if (length(legends) > 0 && length(not_in_data) > 0) {
    warning(sprintf("legends input not in the data: %s", paste(not_in_data, collapse = ", ")))
  }
  legends_in_data <- if (length(legends) == 0) NULL else intersect(legends, present_lineages)
  
  plot_rect <- plot_df %>%
    dplyr::arrange(date, lineage) %>%
    dplyr::group_by(date) %>%
    dplyr::mutate(
      cs_raw    = cumsum(share),
      scale_fac = dplyr::last(cs_raw),
      cs_scaled = dplyr::if_else(scale_fac > 0, cs_raw / scale_fac, 0),
      cs_scaled = pmin(pmax(cs_scaled, 0), 1),
      ymin = dplyr::lag(cs_scaled, default = 0),
      ymax = cs_scaled
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(xmin = date - 0.45,
                  xmax = date + 0.45)
  
  start_date <- as.Date(start_date); end_date <- as.Date(end_date)
  plot_rect  <- plot_rect %>% dplyr::filter(date >= start_date, date <= end_date)
  
  legend_breaks <- if (is.null(legends_in_data) || length(legends_in_data) == 0) {
    NULL
  } else legends_in_data
  
  ggplot2::ggplot(plot_rect) +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = lineage),
      colour = NA
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::scale_x_date(
      limits = c(start_date, end_date),
      expand = c(0, 0),
      breaks = function(x) {
        year_breaks <- seq(as.Date(paste0(format(min(x), "%Y"), "-01-01")),
                           as.Date(paste0(format(max(x), "%Y"), "-01-01")), by = "1 year")
        sort(unique(c(year_breaks, start_date, end_date)))
      },
      labels = function(dates) ifelse(dates %in% c(start_date, end_date),
                                      format(dates, "%Y-%m-%d"),
                                      format(dates, "%Y"))
    ) +
    ggplot2::scale_fill_manual(
      values = cols,
      breaks = legend_breaks,
      guide  = if (is.null(legend_breaks)) "none" else ggplot2::guide_legend(ncol = legend_ncol, title = legend_title)
    ) +
    ggplot2::labs(title = cty, x = NULL, y = "Domestic share") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(linewidth = 0.4, color = "black"),
      axis.ticks       = ggplot2::element_line(linewidth = 0.3, color = "black"),
      axis.ticks.length= grid::unit(2, "mm"),
      axis.text.x      = ggplot2::element_text(size = 9, color = "black", hjust = 0.5),
      axis.text.y      = ggplot2::element_text(size = 9, color = "black"),
      axis.title.y     = ggplot2::element_text(size = 10, color = "black",
                                               margin = ggplot2::margin(r = 6)),
      plot.title       = ggplot2::element_text(size = 11, color = "black", hjust = 0, vjust = 1),
      legend.position  = if (is.null(legend_breaks)) "none" else "right",
      legend.title     = ggplot2::element_text(size = 10, color = "black"),
      legend.text      = ggplot2::element_text(size = 9, color = "black"),
      legend.key.height= grid::unit(4, "mm"),
      legend.key.width = grid::unit(4, "mm"),
      legend.box.margin= ggplot2::margin(0, 0, 0, 6),
      plot.margin      = ggplot2::margin(6, 6, 6, 6)
    )
}

build_variant_panels <- function(data,
                                 countries   = c("United States", "United Kingdom", "South Korea"),
                                 k           = 7,
                                 start_date  = as.Date("2020-01-01"),
                                 end_date    = as.Date("2024-12-31"),
                                 legends_list = list(),
                                 legend_title = NULL,
                                 seed        = NA_integer_) {
  
  make_sample_data <- function(dat, cty) {
    dat %>%
      dplyr::filter(.data$country == cty) %>%
      dplyr::select(country, lineage, date, numerator, denominator, share) %>%
      dplyr::mutate(date = as.Date(date)) %>%
      dplyr::filter(.data$lineage != "Unassigned")
  }
  
  sample_list <- lapply(countries, function(cty) make_sample_data(data, cty))
  names(sample_list) <- countries
  
  all_lineages <- sort(unique(unlist(lapply(sample_list, function(df) unique(df$lineage)))))
  if (!is.na(seed)) set.seed(seed)
  cols_all <- randomcoloR::distinctColorPalette(length(all_lineages))
  names(cols_all) <- all_lineages
  cols_all["others"] <- "lightgray"
  
  cols_list <- lapply(sample_list, function(df) {
    lins <- sort(unique(df$lineage))
    out  <- cols_all[names(cols_all) %in% lins]
    out["others"] <- "lightgray"
    out
  })
  
  plot_list <- mapply(function(df, cols, cty) {
    legends_cty <- legends_list[[cty]]
    plot_variant_stack_smoothed(
      sample_data  = df,
      cols         = cols,
      legend_ncol  = 1,
      k            = k,
      start_date   = start_date,
      end_date     = end_date,
      legends      = if (is.null(legends_cty)) character(0) else legends_cty,
      legend_title = legend_title
    ) + ggplot2::ggtitle(cty)
  }, sample_list, cols_list, names(sample_list), SIMPLIFY = FALSE)
  
  p_all <- Reduce(`/`, plot_list)
  
  list(
    plot_combined     = p_all,
    plots_by_country  = plot_list,
    colors_global     = cols_all,
    colors_by_country = cols_list
  )
}


group_countries_by_volatility_index <- function(data,
                                                n_groups = 3,
                                                label_sep = "-",
                                                weights = c(
                                                  vol_mean_abs      = 0.30,
                                                  vol_sd            = 0.25,
                                                  vol_iqr           = 0.20,
                                                  tail_range        = 0.15,
                                                  vol_entropy       = 0.10,
                                                  pos_neg_imbalance = 0.00
                                                ),
                                                breaks_method = c("jenks", "hclust"),
                                                seed = 123) {
  breaks_method <- match.arg(breaks_method)
  
  req <- c("country_label", "date", "share")
  miss <- setdiff(req, names(data))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  if (!inherits(data$date, "Date")) data$date <- as.Date(data$date)
  
  df <- tidyr::separate(
    data, country_label,
    into = c("country", "classified_label"),
    sep = label_sep, extra = "merge", fill = "right", remove = TRUE
  )
  
  dshare_tbl <- df %>%
    dplyr::arrange(country, classified_label, date) %>%
    dplyr::group_by(country, classified_label) %>%
    dplyr::mutate(dshare = share - dplyr::lag(share)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(is.finite(dshare))
  if (nrow(dshare_tbl) == 0) stop("No finite dshare values found after differencing.")
  
  .safe_entropy <- function(x, nbins = 12L) {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(0)
    rng <- range(x, na.rm = TRUE)
    if (!all(is.finite(rng)) || rng[1] == rng[2]) return(0)
    h <- tryCatch(graphics::hist(x, breaks = nbins, plot = FALSE), error = function(e) NULL)
    if (is.null(h)) return(0)
    p <- h$counts / sum(h$counts); p <- p[p > 0]
    -sum(p * log(p))
  }
  
  country_feat <- dshare_tbl %>%
    dplyr::group_by(country) %>%
    dplyr::summarise(
      n_obs        = dplyr::n(),
      vol_median   = stats::median(dshare, na.rm = TRUE),
      vol_iqr      = stats::IQR(dshare, na.rm = TRUE),
      vol_p10      = stats::quantile(dshare, 0.10, na.rm = TRUE, names = FALSE, type = 7),
      vol_p90      = stats::quantile(dshare, 0.90, na.rm = TRUE, names = FALSE, type = 7),
      vol_mean_abs = mean(abs(dshare), na.rm = TRUE),
      vol_sd       = stats::sd(dshare, na.rm = TRUE),
      pos_ratio    = mean(dshare >  0, na.rm = TRUE),
      zero_ratio   = mean(dshare == 0, na.rm = TRUE),
      neg_ratio    = mean(dshare <  0, na.rm = TRUE),
      vol_entropy  = .safe_entropy(dshare, nbins = 12L),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      tail_range        = vol_p90 - vol_p10,
      pos_neg_imbalance = abs(pos_ratio - neg_ratio)
    )
  
  all_countries <- df %>% dplyr::distinct(country)
  country_feat  <- dplyr::right_join(country_feat, all_countries, by = "country")
  
  num_cols <- setdiff(names(country_feat), "country")
  country_feat[num_cols] <- lapply(country_feat[num_cols], function(x) {
    x[!is.finite(x)] <- NA_real_; x[is.na(x)] <- 0; x
  })
  
  rob_scale <- function(x) {
    med <- stats::median(x, na.rm = TRUE)
    mad <- stats::mad(x, center = med, constant = 1, na.rm = TRUE)
    if (!is.finite(mad) || mad == 0) {
      sdv <- stats::sd(x, na.rm = TRUE)
      if (!is.finite(sdv) || sdv == 0) return(x - med)
      return((x - med) / sdv)
    }
    (x - med) / mad
  }
  
  features_for_index <- intersect(
    c("vol_mean_abs","vol_sd","vol_iqr","tail_range","vol_entropy","pos_neg_imbalance"),
    names(country_feat)
  )
  weights <- weights[features_for_index]
  if (any(!is.finite(weights))) stop("Invalid weights supplied.")
  if (sum(weights) != 0) weights <- weights / sum(weights)
  
  scaled <- country_feat[features_for_index]
  scaled[] <- lapply(scaled, rob_scale)
  scaled[] <- lapply(scaled, function(x) { x[!is.finite(x)] <- 0; x })
  
  volatility_index <- as.numeric(as.matrix(scaled) %*% as.numeric(weights))
  if (!is.finite(sum(volatility_index))) {
    set.seed(seed); volatility_index <- volatility_index + rnorm(length(volatility_index), sd = 1e-6)
  }
  
  tmp <- dplyr::bind_cols(country_feat["country"], scaled) %>%
    dplyr::mutate(volatility_index = volatility_index)
  
  if (breaks_method == "jenks" && requireNamespace("classInt", quietly = TRUE)) {
    ci <- classInt::classIntervals(tmp$volatility_index, n = n_groups, style = "jenks")
    brks <- unique(ci$brks)
    if (length(brks) - 1L < n_groups) {
      set.seed(seed)
      xj <- jitter(tmp$volatility_index, factor = 1e-8)
      ci <- classInt::classIntervals(xj, n = n_groups, style = "jenks")
      brks <- unique(ci$brks)
    }
    cut_lbl <- cut(tmp$volatility_index, breaks = brks, include.lowest = TRUE, right = TRUE)
    grp <- as.integer(factor(cut_lbl, levels = levels(cut_lbl), ordered = TRUE))
    grp_factor <- factor(grp, levels = seq_len(max(grp, na.rm = TRUE)),
                         labels = as.character(seq_len(max(grp, na.rm = TRUE))))
  } else {
    set.seed(seed)
    ord <- order(tmp$volatility_index)
    x   <- tmp$volatility_index[ord]
    hc  <- stats::hclust(stats::dist(x), method = "complete")
    k   <- min(n_groups, length(unique(x)))
    grp_ord <- stats::cutree(hc, k = k)
    grp <- integer(length(x)); grp[seq_along(x)] <- grp_ord
    grp <- grp[order(order(tmp$volatility_index))]  # align to tmp order
    grp_factor <- factor(grp, levels = sort(unique(grp)),
                         labels = as.character(sort(unique(grp))))
  }
  
  out <- dplyr::bind_cols(
    tmp,
    setNames(scaled, paste0(names(scaled), "_z"))
  ) %>%
    dplyr::mutate(volatility_group = grp_factor)
  
  out
}







build_lineage_periods <- function(data2,
                                  WIN     = 14L,
                                  SUM_MIN = 10L,
                                  NZ_MIN  = 3L,
                                  show_progress = TRUE) {
  req <- c("country", "lineage", "date", "numerator")
  miss <- setdiff(req, names(data2))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  if (!inherits(data2$date, "Date")) data2$date <- as.Date(data2$date)
  
  combo_list <- data2 %>% dplyr::distinct(country, lineage)
  n_total <- nrow(combo_list)
  pb <- if (isTRUE(show_progress)) utils::txtProgressBar(min = 0, max = n_total, style = 3) else NULL
  lineage_periods <- vector("list", n_total)
  
  for (i in seq_len(n_total)) {
    cty <- combo_list$country[i]
    lin <- combo_list$lineage[i]
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
    
    df_sub <- data2 %>%
      dplyr::filter(.data$country == cty, .data$lineage == lin) %>%
      dplyr::arrange(.data$date)
    
    n_obs_days      <- dplyr::n_distinct(df_sub$date)
    total_numerator <- sum(df_sub$numerator, na.rm = TRUE)
    
    out_row <- data.frame(
      country         = cty,
      lineage         = lin,
      start_date      = as.Date(NA),
      end_date        = as.Date(NA),
      n_obs_days      = n_obs_days,
      total_numerator = total_numerator,
      selected        = FALSE,
      reason          = NA_character_,
      stringsAsFactors = FALSE
    )
    
    if (nrow(df_sub) < WIN) {
      out_row$start_date <- min(df_sub$date, na.rm = TRUE)
      out_row$end_date   <- max(df_sub$date, na.rm = TRUE)
      out_row$reason     <- sprintf("too_short: fewer than %d daily observations", WIN)
      lineage_periods[[i]] <- out_row
      next
    }
    
    df_sub <- df_sub %>%
      dplyr::mutate(
        roll_sum = zoo::rollapply(.data$numerator, width = WIN, FUN = sum, fill = NA, align = "left"),
        roll_nonzero_days = zoo::rollapply(.data$numerator > 0, width = WIN, FUN = sum, fill = NA, align = "left")
      )
    
    valid_rows <- df_sub %>%
      dplyr::filter(.data$roll_sum >= SUM_MIN, .data$roll_nonzero_days >= NZ_MIN) %>%
      dplyr::select(date)
    
    if (nrow(valid_rows) == 0) {
      cond_sum_ok <- any(df_sub$roll_sum >= SUM_MIN, na.rm = TRUE)
      cond_nz_ok  <- any(df_sub$roll_nonzero_days >= NZ_MIN, na.rm = TRUE)
      if (!cond_sum_ok && !cond_nz_ok) {
        out_row$reason <- sprintf("no_window: sum<%d AND nonzero_days<%d", SUM_MIN, NZ_MIN)
      } else if (!cond_sum_ok) {
        out_row$reason <- sprintf("no_window: sum<%d", SUM_MIN)
      } else if (!cond_nz_ok) {
        out_row$reason <- sprintf("no_window: nonzero_days<%d", NZ_MIN)
      } else {
        out_row$reason <- "no_window: mixed_failure"
      }
      lineage_periods[[i]] <- out_row
      next
    }
    
    dseq <- sort(as.Date(valid_rows$date))
    grp  <- c(1L, 1L + cumsum(diff(dseq) > 1))
    runs <- data.frame(
      grp        = grp,
      date       = dseq,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::group_by(grp) %>%
      dplyr::summarise(
        run_start = min(date),
        run_end   = max(date),
        run_len   = dplyr::n(),
        .groups   = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(run_len), dplyr::desc(run_end))
    
    chosen <- runs[1, , drop = FALSE]
    
    out_row$start_date <- chosen$run_start
    out_row$end_date   <- chosen$run_end + (WIN - 1L)
    out_row$selected   <- TRUE
    out_row$reason     <- "passed"
    
    lineage_periods[[i]] <- out_row
  }
  
  if (!is.null(pb)) close(pb)
  lineage_period_df <- dplyr::bind_rows(lineage_periods)
  tibble::as_tibble(lineage_period_df)
}






plot_top_lineages_with_periods <- function(country_name,
                                           df,
                                           analytical_period_df,
                                           n = 9,
                                           start_date = NULL,
                                           end_date   = NULL) {
  stopifnot(all(c("country","lineage","date","share") %in% names(df)))
  stopifnot(all(c("country","analytical_start","analytical_end") %in% names(analytical_period_df)))
  
  if (!inherits(df$date, "Date")) df$date <- as.Date(df$date)
  if (!inherits(analytical_period_df$analytical_start, "Date")) {
    analytical_period_df$analytical_start <- as.Date(analytical_period_df$analytical_start)
  }
  if (!inherits(analytical_period_df$analytical_end, "Date")) {
    analytical_period_df$analytical_end <- as.Date(analytical_period_df$analytical_end)
  }
  
  df_cty <- df[df$country == country_name, , drop = FALSE]
  if (nrow(df_cty) == 0) stop("No rows for country: ", country_name)
  
  ap <- analytical_period_df[analytical_period_df$country == country_name, , drop = FALSE]
  if (nrow(ap) == 0) {
    warning("No analytical period found for country: ", country_name, ". Lines will be omitted.")
    ap <- NULL
  }
  
  top_lin <- df_cty |>
    dplyr::group_by(lineage) |>
    dplyr::summarise(sum_share = sum(share, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(sum_share)) |>
    dplyr::slice_head(n = n) |>
    dplyr::pull(lineage)
  
  # special lineages you want fixed at rightmost of row 1 and 2
  special <- c("JN.1.4", "JN.1.16")
  available_special <- intersect(special, top_lin)
  others <- setdiff(top_lin, available_special)
  
  # 构造 grid 顺序
  ncol <- 3L
  slots <- vector("list", length = length(top_lin))
  
  row1_idx <- seq(1, ncol)       # 第一行位置
  row2_idx <- seq(ncol + 1, 2*ncol) # 第二行位置
  
  if ("JN.1.4" %in% available_special) {
    slots[[row1_idx[ncol]]] <- "JN.1.4"
    available_special <- setdiff(available_special, "JN.1.4")
  }
  if ("JN.1.16" %in% available_special) {
    slots[[row2_idx[ncol]]] <- "JN.1.16"
    available_special <- setdiff(available_special, "JN.1.16")
  }
  
  # 填剩下的空位
  fill_order <- c(others, available_special)
  empty_positions <- which(sapply(slots, is.null))
  slots[empty_positions] <- fill_order[seq_along(empty_positions)]
  
  top_lin <- unlist(slots)
  
  n_panels <- length(top_lin)
  nrow <- ceiling(n_panels / ncol)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(nrow, ncol), mar = c(3, 3, 2, 1))
  
  for (lin in top_lin) {
    dl <- df_cty[df_cty$lineage == lin, c("date","share"), drop = FALSE]
    dl <- dl[order(dl$date), , drop = FALSE]
    
    if (nrow(dl) == 0) {
      plot.new(); title(main = lin); next
    }
    
    xlim <- NULL
    if (!is.null(start_date) && !is.null(end_date)) {
      xlim <- c(as.Date(start_date), as.Date(end_date))
    }
    
    plot(dl$date, dl$share,
         type = "p", ylim = c(0, 1),
         xlab = "", ylab = "Share",
         main = lin,
         xlim = xlim,
         pch = 3,
         cex = 0.6)
    
    if (!is.null(ap)) {
      abline(v = ap$analytical_start[1], col = "orange", lty = 1, lwd = 2)
      abline(v = ap$analytical_end[1],   col = "orange", lty = 1, lwd = 2)
      
      text(x = ap$analytical_start[1], y = 0.9,
           labels = format(ap$analytical_start[1], "%b %Y"),
           pos = 4, col = "orange", font = 2, cex = 0.8)
      
      text(x = ap$analytical_end[1], y = 0.9,
           labels = format(ap$analytical_end[1], "%b %Y"),
           pos = 2, col = "orange", font = 2, cex = 0.8)
    }
  }
  
  invisible(TRUE)
}





split_one <- function(feat_df, meas_df) {
  df <- feat_df %>%
    dplyr::inner_join(meas_df, by = c("country","lineage"))
  
  id_cols   <- c("country","lineage")
  y_col     <- "peak_share"
  feat_cols <- setdiff(names(df), c(id_cols, y_col))
  
  df[feat_cols] <- lapply(df[feat_cols], as.numeric)
  
  all_na_row <- apply(df[feat_cols], 1, function(r) all(is.na(r)))
  if (any(all_na_row)) df <- df[!all_na_row, , drop = FALSE]
  
  X <- as.matrix(df[, feat_cols, drop = FALSE])
  Y <- as.numeric(df[[y_col]])
  
  list(
    X = X,
    Y = Y,
    ids = df[, id_cols, drop = FALSE],
    feature_names = feat_cols
  )
}
