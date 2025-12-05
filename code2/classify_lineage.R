


# Function: get_classification
# Purpose: Iteratively group related SARS-CoV-2 lineages based on frequency and count thresholds,
#          then map them to alias names for simplified classification.
# Inputs:
#   - lineages: vector of lineage names (e.g., "BA.2.75.2")
#   - number_sequences: vector of counts for each lineage
#   - p_lim: frequency threshold (default = 0.01)
#   - n_lim: count threshold (default = 50)
#   - alias_list: alias mapping table from get_alias()
# Output:
#   - tibble containing final classified results with alias and label columns

get_classification <- function(lineages, number_sequences, p_lim = 0.01, alias_list) {
  
  # Initial lineage table with proportions and flags for threshold filtering
  classifications <- tidyr::tibble(
    lineages = lineages,
    n = number_sequences,
    N = sum(number_sequences)
  ) |>
    dplyr::mutate(
      p = n / N
    ) |>
    dplyr::mutate(
      classified_interim = lineages,
      classified_final = dplyr::if_else(p > p_lim,
                                        lineages,
                                        NA)
    )
  
  dim_old <- nrow(classifications)
  
  # Iteratively trace back to parent lineages until thresholds are met
  while (any(is.na(classifications$classified_final))) {
    
    classifications |>
      dplyr::pull(classified_interim) -> possible_parents
    
    # Identify possible parent-child structure based on lineage string prefixes
    found_in <- lapply(possible_parents, stringr::str_starts, pattern = paste0(possible_parents, "."))
    
    classifications |>
      dplyr::rowwise() |>
      dplyr::mutate(
        # Find all higher-level parent lineages that current lineage starts with
        best_matches = list(possible_parents[stringr::str_starts(classified_interim,
                                                                 stringr::coll(paste0(possible_parents, ".")))])
      ) |>
      dplyr::mutate(
        # Choose the longest matching parent lineage (most specific)
        selected_parent = ifelse(length(best_matches) == 0,
                                 NA,
                                 best_matches[which.max(stringr::str_length(best_matches))]),
        # Remove the last segment if no matching parent found (go one level up)
        cropped = stringr::str_split(classified_interim, stringr::coll(".")),
        cropped_processed = ifelse(length(cropped) > 2,
                                   paste((cropped)[-length((cropped))], collapse = "."),
                                   paste((cropped), collapse = ".")),
        # Update lineage label for next iteration
        new_list = ifelse(is.na(selected_parent),
                          cropped_processed,
                          selected_parent)
      )  |>
      dplyr::ungroup() |>
      dplyr::mutate(
        # Replace NA classifications with the new upper-level lineage
        classified_interim = dplyr::if_else(is.na(classified_final),
                                            new_list,
                                            classified_final)
      ) |>
      # Re-aggregate by current lineage group
      dplyr::summarise(
        n = sum(n),
        .by = c("classified_interim", "N")
      ) |>
      dplyr::mutate(
        p = n / N
      ) |>
      # Update classification if threshold is met (only p_lim)
      dplyr::mutate(
        classified_final = dplyr::if_else(p > p_lim,
                                          classified_interim,
                                          NA)
      ) -> classifications
    
    # Stop if no further aggregation occurs
    if (nrow(classifications) == dim_old) {
      break("Hunting over")
    }
    
    dim_old <- nrow(classifications)
  }
  
  # Map final classifications to alias names
  classifications |>
    tidyr::drop_na() |>
    fuzzyjoin::fuzzy_left_join(
      alias_list,
      by = c("classified_final" = "alias_lineage"),
      match_fun = stringr::str_starts
    ) |>
    dplyr::mutate(
      alias_lineage = dplyr::if_else(is.na(alias_lineage),
                                     classified_final,
                                     alias_lineage),
      length_alias = stringr::str_length(alias_lineage)
    ) |>
    # Keep the longest alias match for each lineage
    dplyr::filter(
      length_alias == max(length_alias),
      .by = classified_final
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      # Construct alias-based naming convention
      lineage_class_alias = gsub(x = classified_final,
                                 pattern = alias_lineage,
                                 replacement = "",
                                 fixed = TRUE),
      classified_alias = dplyr::if_else(is.na(alias),
                                        classified_final,
                                        paste0(alias, lineage_class_alias)),
      classified_label = dplyr::if_else(is.na(alias),
                                        classified_final,
                                        paste0(classified_alias, " (", classified_final, ")"))
    ) |>
    dplyr::rename(
      classified_unasliased = classified_final
    ) |>
    dplyr::distinct(
      classified_unasliased,
      classified_alias,
      classified_label
    ) |>
    dplyr::ungroup() -> reference_frame
  
  return(reference_frame)
}





classify_and_summarise_variants <- function(data_cty, alias_list, p_lim = .01) {
  
  # p_lim = 0.01; alias_list = alias_map; # for test
  # data_cty <- data3 %>%
  #   filter(country == "Japan",
  #          !is.na(lineage),
  #          lineage != "")
  
  global_wave_periods <- list(
    wave1 = c("2020-03-01", "2020-06-30"),  # 原始株第一波
    wave2 = c("2020-10-01", "2021-02-28"),  # 秋冬第二波 / Alpha 前后
    wave3 = c("2021-03-01", "2021-11-30"),  # Delta 波
    wave4 = c("2021-12-01", "2022-03-31"),  # Omicron BA.1/BA.2
    wave5 = c("2022-05-01", "2022-09-30"),  # Omicron BA.4/BA.5 夏波
    wave6 = c("2022-11-01", "2024-04-30")   # XBB/BQ 等后 Omicron 波
  )
  
  summarised <- lapply(global_wave_periods, function(x) {
    df <- data_cty %>%
      filter(date >= as.Date(x[1]), date <= as.Date(x[2])) %>%
      group_by(lineage) %>%
      summarise(n = sum(numerator), .groups = "drop")
    if (nrow(df) > 0) df else NULL
  })
  
  names(summarised) <- names(global_wave_periods)
  summarised <- summarised %>% purrr::compact() 
  
  
  classification_to_append <- lapply(names(summarised), function(w) {
    get_classification(
      lineages = summarised[[w]]$lineage,
      number_sequences = summarised[[w]]$n,
      p_lim = p_lim,
      alias_list = alias_list
    ) %>%
      mutate(wave = w,
             decimal_lineage = paste0(classified_unasliased, "."))
  }) %>%
    bind_rows() %>% 
    distinct(classified_unasliased, classified_alias, classified_label, .keep_all = TRUE)
  
  data_cty_classfied <- data_cty %>%
    dplyr::mutate(decimal_lineage = paste0(lineage, ".")) %>%
    fuzzyjoin::fuzzy_left_join(
      classification_to_append,
      by = "decimal_lineage",
      match_fun = stringr::str_starts
    ) %>%
    dplyr::filter(
      is.na(classified_unasliased) |
        stringr::str_length(classified_unasliased) ==
        max(stringr::str_length(classified_unasliased)),
      .by = lineage
    ) %>%
    dplyr::distinct(lineage, date, .keep_all = TRUE) %>%
    dplyr::mutate(
      length_class = stringr::str_count(classified_unasliased, "\\."),
      length_lineage = stringr::str_count(lineage, "\\."),
      classified_label = dplyr::if_else(
        (length_lineage - length_class) <= 1,
        classified_label,
        NA_character_
      ),
      classified_label = dplyr::if_else(
        is.na(classified_label), "Other", classified_label
      )
    ) %>%
    dplyr::filter(classified_label != "Other") %>%
    dplyr::select(country, date, lineage, numerator, denominator, share, classified_label)
  
  data_cty_classfied <- data_cty_classfied %>%
    dplyr::group_by(date, classified_label) %>%
    dplyr::summarise(
      numerator = sum(numerator, na.rm = TRUE),
      share = sum(share, na.rm = TRUE), 
      denominator = dplyr::first(denominator),
      .groups = "drop"
    )
  
  list(
    data = data_cty_classfied,
    classification = classification_to_append
  )
}



# # load the data 
# data3 <- readRDS("code/data3_and_dropped.rds")[["data3"]]

process_classified_variants <- function(data, alias_list,
                                        p_lim = 0.01) {
  
  # data = data3; p_lim = 0.01; alias_list = alias_map # for test
  
  country_list <- unique(data$country)
  classified_result_list <- list()

  for (cty in country_list) {
    
    # cty = "South Korea" # for test
    
    message("Processing: ", cty)
    
    data_cty <- data %>%
      filter(country == cty,
             !is.na(lineage),
             lineage != "")
    
    classified_result_list[[cty]] <- classify_and_summarise_variants(
      data_cty, alias_list, p_lim = p_lim
    )
  }

  nested_classified_list <- list()
  for (cty in country_list) {
    df <- classified_result_list[[cty]]$data
    denom_by_date <- df %>% dplyr::distinct(date, denominator)
    
    nested_classified_list[[cty]] <- list()
    label_list <- unique(df$classified_label)
    
    for (label in label_list) {
      df_label <- df %>%
        dplyr::filter(classified_label == label) %>%
        dplyr::arrange(date)
      
      if (nrow(df_label) == 0) next
      
      date_range <- range(df_label$date, na.rm = TRUE)
      all_dates <- seq.Date(from = date_range[1], to = date_range[2], by = "day")
      
      complete_df <- tidyr::expand_grid(date = all_dates, classified_label = label) %>%
        dplyr::left_join(denom_by_date, by = "date") %>%
        dplyr::left_join(df_label, by = c("date", "classified_label")) %>%
        dplyr::mutate(
          country = cty,
          numerator = ifelse(is.na(numerator), 0, numerator),
          share = ifelse(is.na(share), 0, share),
          denominator = dplyr::coalesce(denominator.y, denominator.x)
        ) %>%
        dplyr::select(country, date, classified_label, numerator, denominator, share) %>%
        dplyr::arrange(date)
      
      nested_classified_list[[cty]][[label]] <- complete_df
    }
  }
  
  # Step 3: flatten nested list into one big tibble
  flat_df_list <- list()
  for (cty in names(nested_classified_list)) {
    for (label in names(nested_classified_list[[cty]])) {
      df <- nested_classified_list[[cty]][[label]]
      flat_df_list[[paste(cty, label, sep = "|")]] <- df
    }
  }
  
  classified_data <- dplyr::bind_rows(flat_df_list) %>%
    dplyr::mutate(country_label = paste0(country, "-", classified_label)) %>%
    dplyr::select(country_label, date, numerator, denominator, share)
  
  return(classified_data)
}


