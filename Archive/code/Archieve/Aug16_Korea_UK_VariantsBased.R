# source the functions (functions_0706.R)
source('Background/functions_0706.R')
source('Code/functions_franky_0725.R')
source('Background/growth_rate_function_original.R')

# load packages 
library(magrittr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(viridis)
library(scales) 
library(fda)
library(tidyr)

# get the alias
alias_list <- get_alias()

# set the time range (as large as possible)
months_from_most_recent <- 365 * 6

# load the data and process 
url <- 'Data/summary_gisaid_new.csv'
data <- readr::read_csv(url) |>
  dplyr::filter(
    # country == country_select,
    lineage != 'Unassigned',
    date >= max(date, na.rm = T) - months_from_most_recent
  ) |>
  dplyr::mutate(
    alias = gsub(
      pattern = "[^a-zA-Z]",
      replacement = "",
      x = lineage
    ),
    second_part = gsub(pattern = "[a-zA-Z]",
                       replacement = "",
                       x = lineage)
  ) |>
  dplyr::left_join(alias_list) |>
  dplyr::mutate(
    unaliased_lineage  =
      dplyr::if_else(!is.na(alias_lineage),
                     paste0(alias_lineage, second_part),
                     lineage),
  )

dataset <- data %>%
  dplyr::select(country, date, lineage, numerator, unaliased_lineage) %>%
  dplyr::rename(n = numerator)

# classified the variants
summarised <- dataset %>%
  dplyr::group_by(unaliased_lineage) %>%
  dplyr::summarise(n = sum(n))

# # set p limit
# p_lim <- 0.001
# 
# # set n limit
# n_lim <- 100
# # get the classification to append
# # save the result since the function costs a lot of time
# classification_to_append <- get_classification(
#   lineages = summarised$unaliased_lineage,
#   number_sequences = summarised$n,
#   p_lim = p_lim, n_lim = n_lim, alias_list
# ) |>
#   dplyr::mutate(
#     decimal_lineage = paste0(classified_unasliased, '.')
#   )
# 
# saveRDS(classification_to_append, "Data/Aug17_classification_to_append_Variants.rds")


# load the .rds file
append_url <- "Data/Aug17_classification_to_append_Variants.rds"
classification_to_append <- readRDS(append_url)

# append the classification
classified_data <- dataset |>
  dplyr::mutate(
    decimal_lineage = paste0(unaliased_lineage, '.')
  ) |>
  fuzzyjoin::fuzzy_left_join(
    classification_to_append,
    by = 'decimal_lineage',
    match_fun = stringr::str_starts
  ) |>
  dplyr::filter(
    is.na(classified_unasliased) | stringr::str_length(classified_unasliased) == max(stringr::str_length(classified_unasliased)),
    .by = c(lineage, date)
  ) |>
  dplyr::mutate(
    length_class = stringr::str_count(classified_unasliased, '\\.'),
    length_lineage = stringr::str_count(unaliased_lineage, '\\.'),
    classified_label = dplyr::if_else((length_lineage - length_class) <= 1,
                                      classified_label,
                                      NA),
    classified_label = dplyr::if_else(is.na(classified_label),
                                      "Other",
                                      classified_label)
  ) |>
  dplyr::filter(classified_label != 'Other') |>
  dplyr::select(date, n, classified_label)

# # extract the group
# classified_data$Group <- ifelse(grepl("XBB", classified_data$classified_label), "XBB",
#                                 ifelse(grepl("BA", classified_data$classified_label), "BA", "B"))

classified_data %>%
  complete_dates() %>%
  dplyr::mutate(
    Group = ifelse(grepl("XBB", classified_label), "XBB",
                   ifelse(grepl("BA", classified_label), "BA", "B"))
  ) %>%
  dplyr::group_by(date, Group, classified_label) %>%
  dplyr::summarise(numerator = sum(n), .groups = 'drop') %>%
  dplyr::mutate(
    denominator = sum(numerator, na.rm = TRUE),
    share = numerator / denominator,
    share = ifelse(is.na(share), 0, share)
  ) %>%
  dplyr::group_by(classified_label) %>%
  dplyr::mutate(
    gap = as.numeric(date - dplyr::first(date))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Group, classified_label) -> variants_data

variants_data_split <- split(variants_data, variants_data$Group)

variants_data %>% 
  dplyr::group_by(
    classified_label
  ) %>% 
  dplyr::summarise(
    response = sum(share), 
    .groups = 'drop'
  ) %>% 
  arrange(
    desc(response)
  ) -> variants_dominance


rm(list = setdiff(ls(), c("variants_data_split", "classified_data", "variants_data", "variants_dominance",
                          "extract_growth_rates", "paa_transform")))

variants_data %>% 
  dplyr::group_by(
    Group, classified_label
  ) %>% 
  dplyr::summarise(
    n_total <- sum(numerator)
  ) -> summary_table


generate_growth_rate_list <- function(df) {
  
  # Function to process growth rates for a given column
  process_growth_rates <- function(data, column_name) {
    output_ <- extract_growth_rates(data[[column_name]], data$date, dow = "none", denominator = 3)
    list(
      fit = output_$model_fit %>% 
        dplyr::select(Date_numeric, date, full_fit) %>% 
        dplyr::rename(!!paste0(column_name, "_fit") := full_fit),
      gr_fit = output_$fitted_rates %>% 
        dplyr::select(Date_numeric, date, rate) %>% 
        dplyr::rename(!!paste0(column_name, "_gr_fit") := rate)
    )
  }
  
  # Initialize the output list
  list_ <- list()
  
  # Loop over each variant and process data
  for (variant_ in unique(df$classified_label)) {
    
    cat("processing-", variant_, "\n")
    
    variant_df_ <- df %>% 
      dplyr::filter(classified_label == variant_)
    
    # Process 'numerator' and 'share' columns
    abs_data <- process_growth_rates(variant_df_, "numerator")
    share_data <- process_growth_rates(variant_df_, "share")
    
    # Combine the results
    output_ <- abs_data$fit %>% 
      dplyr::left_join(abs_data$gr_fit, by = c("Date_numeric", "date")) %>% 
      dplyr::left_join(share_data$fit, by = c("Date_numeric", "date")) %>% 
      dplyr::left_join(share_data$gr_fit, by = c("Date_numeric", "date"))
    
    # Store the raw data and the fitted output in the list
    list_[[variant_]] <- list(
      raw = variant_df_,
      fit = output_
    )
  }
  
  # Return the final list
  return(list_)
}



variants_ls_30 <- list()
variants_ls_60 <- list()

for (name in names(variants_data_split)) {
  
  # name <- "B"
  
  cat("currently processing: ", name, "\n")
  
  df_1_ <- variants_data_split[[name]] %>% 
    dplyr::filter(
      gap <= 30
    )
  
  df_2_ <- variants_data_split[[name]] %>% 
    dplyr::filter(
      gap <= 60
    )
  
  variants_ls_30[[name]] <- generate_growth_rate_list(df_1_)
  variants_ls_60[[name]] <- generate_growth_rate_list(df_2_)
  
}

################################################################################










abs_raw_ls_30 <- list()
share_raw_ls_30 <- list()
abs_fit_ls_30 <- list()
abs_gr_ls_30 <- list()
share_fit_ls_30 <- list()
share_gr_ls_30 <- list()

abs_raw_df_30 <- list()
share_raw_df_30 <- list()
abs_fit_df_30 <- list()
abs_gr_df_30 <- list()
share_fit_df_30 <- list()
share_gr_df_30 <- list()


for (key in names(variants_ls_30)) {
  
  # key = "B"
  
  cat("processing- ", key, "\n")
  
  ls_ <- variants_ls_30[[key]]
  
  for (name in names(ls_)) {
    
    # name = "AD.2 (B.1.1.315.2)"
    
    cat("    processing- ", key, ": ", name, "\n")
    
    select_step <- seq(from = 1, to = 301, by = 10)
    
    abs_raw_col <- c(name, ls_[[name]]$raw$numerator %>% as.vector())
    share_raw_col <- c(name, ls_[[name]]$raw$share %>% as.vector())
    
    abs_fit_col <- c(name, (ls_[[name]]$fit$numerator_fit %>% as.vector())[select_step])
    abs_gr_col <- c(name, (ls_[[name]]$fit$numerator_gr_fit %>% as.vector())[select_step])
    
    share_fit_col <- c(name, (ls_[[name]]$fit$share_fit %>% as.vector())[select_step])
    share_gr_col <- c(name, (ls_[[name]]$fit$share_gr_fit %>% as.vector())[select_step])
    
    abs_raw_ls_30[[key]][[name]] <- abs_raw_col
    share_raw_ls_30[[key]][[name]] <- share_raw_col
    abs_fit_ls_30[[key]][[name]] <- abs_fit_col
    abs_gr_ls_30[[key]][[name]] <- abs_gr_col
    share_fit_ls_30[[key]][[name]] <- share_fit_col
    share_gr_ls_30[[key]][[name]] <- share_gr_col
    
  }
  
  abs_raw_df_30[[key]] <- do.call(rbind, abs_raw_ls_30[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  share_raw_df_30[[key]] <- do.call(rbind, share_raw_ls_30[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  abs_fit_df_30[[key]] <- do.call(rbind, abs_fit_ls_30[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  abs_gr_df_30[[key]] <- do.call(rbind, abs_gr_ls_30[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  share_fit_df_30[[key]] <- do.call(rbind, share_fit_ls_30[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  share_gr_df_30[[key]] <- do.call(rbind, share_gr_ls_30[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  
  col_names <- c("classified_label", paste0("Day", 0:30))
  colnames(abs_raw_df_30[[key]]) <- col_names
  colnames(share_raw_df_30[[key]]) <- col_names
  colnames(abs_fit_df_30[[key]]) <- col_names
  colnames(abs_gr_df_30[[key]]) <- col_names
  colnames(share_fit_df_30[[key]]) <- col_names
  colnames(share_gr_df_30[[key]]) <- col_names
  
  abs_raw_df_30[[key]] <- abs_raw_df_30[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  share_raw_df_30[[key]] <- share_raw_df_30[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  abs_fit_df_30[[key]] <- abs_fit_df_30[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  abs_gr_df_30[[key]] <- abs_gr_df_30[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  share_fit_df_30[[key]] <- share_fit_df_30[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  share_gr_df_30[[key]] <- share_gr_df_30[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  
  rownames(abs_raw_df_30[[key]]) <- NULL
  rownames(share_raw_df_30[[key]]) <- NULL
  rownames(abs_fit_df_30[[key]]) <- NULL
  rownames(abs_gr_df_30[[key]]) <- NULL
  rownames(share_fit_df_30[[key]]) <- NULL
  rownames(share_gr_df_30[[key]]) <- NULL
  
  colnames_ <- c('response', '1-5', '6-10', '11-15', '16-20', '21-25', '26-30')
  
  abs_raw_df_30[[key]] <- abs_raw_df_30[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  share_raw_df_30[[key]] <- share_raw_df_30[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  abs_fit_df_30[[key]] <- abs_fit_df_30[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  abs_gr_df_30[[key]] <- abs_gr_df_30[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  share_fit_df_30[[key]] <- share_fit_df_30[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  share_gr_df_30[[key]] <- share_gr_df_30[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  
}

ls_30 <- list()
ls_30[["raw_abs"]] <- abs_raw_df_30
ls_30[["raw_share"]] <- share_raw_df_30
ls_30[["fit_abs"]] <- abs_fit_df_30
ls_30[["fit_abs_gr"]] <- abs_gr_df_30
ls_30[["fit_share"]] <- share_fit_df_30
ls_30[["fit_share_gr"]] <- share_gr_df_30










abs_raw_ls_60 <- list()
share_raw_ls_60 <- list()
abs_fit_ls_60 <- list()
abs_gr_ls_60 <- list()
share_fit_ls_60 <- list()
share_gr_ls_60 <- list()

abs_raw_df_60 <- list()
share_raw_df_60 <- list()
abs_fit_df_60 <- list()
abs_gr_df_60 <- list()
share_fit_df_60 <- list()
share_gr_df_60 <- list()


for (key in names(variants_ls_60)) {
  
  # key = "B"
  
  cat("processing- ", key, "\n")
  
  ls_ <- variants_ls_60[[key]]
  
  for (name in names(ls_)) {
    
    # name = "AD.2 (B.1.1.315.2)"
    
    cat("    processing- ", key, ": ", name, "\n")
    
    select_step <- seq(from = 1, to = 601, by = 10)
    
    abs_raw_col <- c(name, ls_[[name]]$raw$numerator %>% as.vector())
    share_raw_col <- c(name, ls_[[name]]$raw$share %>% as.vector())
    
    abs_fit_col <- c(name, (ls_[[name]]$fit$numerator_fit %>% as.vector())[select_step])
    abs_gr_col <- c(name, (ls_[[name]]$fit$numerator_gr_fit %>% as.vector())[select_step])
    
    share_fit_col <- c(name, (ls_[[name]]$fit$share_fit %>% as.vector())[select_step])
    share_gr_col <- c(name, (ls_[[name]]$fit$share_gr_fit %>% as.vector())[select_step])
    
    abs_raw_ls_60[[key]][[name]] <- abs_raw_col
    share_raw_ls_60[[key]][[name]] <- share_raw_col
    abs_fit_ls_60[[key]][[name]] <- abs_fit_col
    abs_gr_ls_60[[key]][[name]] <- abs_gr_col
    share_fit_ls_60[[key]][[name]] <- share_fit_col
    share_gr_ls_60[[key]][[name]] <- share_gr_col
    
  }
  
  abs_raw_df_60[[key]] <- do.call(rbind, abs_raw_ls_60[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  share_raw_df_60[[key]] <- do.call(rbind, share_raw_ls_60[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  abs_fit_df_60[[key]] <- do.call(rbind, abs_fit_ls_60[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  abs_gr_df_60[[key]] <- do.call(rbind, abs_gr_ls_60[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  share_fit_df_60[[key]] <- do.call(rbind, share_fit_ls_60[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  share_gr_df_60[[key]] <- do.call(rbind, share_gr_ls_60[[key]]) %>% as.data.frame(stringsAsFactors = FALSE)
  
  col_names <- c("classified_label", paste0("Day", 0:60))
  colnames(abs_raw_df_60[[key]]) <- col_names
  colnames(share_raw_df_60[[key]]) <- col_names
  colnames(abs_fit_df_60[[key]]) <- col_names
  colnames(abs_gr_df_60[[key]]) <- col_names
  colnames(share_fit_df_60[[key]]) <- col_names
  colnames(share_gr_df_60[[key]]) <- col_names
  
  abs_raw_df_60[[key]] <- abs_raw_df_60[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  share_raw_df_60[[key]] <- share_raw_df_60[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  abs_fit_df_60[[key]] <- abs_fit_df_60[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  abs_gr_df_60[[key]] <- abs_gr_df_60[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  share_fit_df_60[[key]] <- share_fit_df_60[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  share_gr_df_60[[key]] <- share_gr_df_60[[key]] %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
  
  rownames(abs_raw_df_60[[key]]) <- NULL
  rownames(share_raw_df_60[[key]]) <- NULL
  rownames(abs_fit_df_60[[key]]) <- NULL
  rownames(abs_gr_df_60[[key]]) <- NULL
  rownames(share_fit_df_60[[key]]) <- NULL
  rownames(share_gr_df_60[[key]]) <- NULL
  
  colnames_ <- c('response', '1-5', '6-10', '11-15', '16-20', '21-25', '26-30', 
                 '31-35', '36-40', '41-45', '46-50', '51-55', '55-60')
  
  # colnames_ <- c('response', '1-5', '6-10', '11-15', '16-20', '21-25', '26-30')
  
  abs_raw_df_60[[key]] <- abs_raw_df_60[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  share_raw_df_60[[key]] <- share_raw_df_60[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  abs_fit_df_60[[key]] <- abs_fit_df_60[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  abs_gr_df_60[[key]] <- abs_gr_df_60[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  share_fit_df_60[[key]] <- share_fit_df_60[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  share_gr_df_60[[key]] <- share_gr_df_60[[key]] %>% dplyr::left_join(variants_dominance, by = "classified_label")
  
}


ls_60 <- list()
ls_60[["raw_abs"]] <- abs_raw_df_60
ls_60[["raw_share"]] <- share_raw_df_60
ls_60[["fit_abs"]] <- abs_fit_df_60
ls_60[["fit_abs_gr"]] <- abs_gr_df_60
ls_60[["fit_share"]] <- share_fit_df_60
ls_60[["fit_share_gr"]] <- share_gr_df_60

################################################################################






ls_30_swapped <- list()

for (level1_key in names(ls_30)) {
  for (level2_key in names(ls_30[[level1_key]])) {
    if (is.null(ls_30_swapped[[level2_key]])) {
      ls_30_swapped[[level2_key]] <- list()
    }
    ls_30_swapped[[level2_key]][[level1_key]] <- ls_30[[level1_key]][[level2_key]]
  }
}


# generate models
model_30days_data <- list()

for (group in names(ls_30_swapped)) {
  
  colnames_ <- c('response', '1-5', '6-10', '11-15', '16-20', '21-25', '26-30')
  
  # group <- "B"
  
  cat("-fitting: ", group, "\n")
  
  for (type in names(ls_30_swapped[[group]])) {
    
    # type = "raw_abs"
    # type = "raw_share"
    
    cat("    -processing: ", type, "\n")
    
    Y <- ls_30_swapped[[group]][[type]]$response
    X <- ls_30_swapped[[group]][[type]][, paste0("Day", 0:30)] %>% 
      as.matrix() %>% 
      t() %>% 
      apply(2, paa_transform, num_segments = 6) %>% 
      t()
      
    input_df <- cbind(Y, X) %>% as.data.frame()
    colnames(input_df) <- colnames_
    
    model_30days_data[[group]][[type]] <- input_df
    
  }
  
}

model_60days <- list()

for (type in names(ls_60)) {
  
  colnames_ <- c('response', '1-5', '6-10', '11-15', '16-20', '21-25', '26-30', 
                 '31-35', '36-40', '41-45', '46-50', '51-55', '55-60')
  
  # type <- "raw_abs"
  
  cat("-fitting: ", type, "\n")
  
  for (group in names(ls_60[[type]])) {
    
    # group <- "BA"
    
    cat("    -processing: ", group, "\n")
    
    Y <- ls_60[[type]][[group]]$response
    X <- ls_60[[type]][[group]][, paste0("Day", 0:60)] %>% 
      as.matrix() %>% 
      t() %>% 
      apply(2, paa_transform, num_segments = 12) %>% 
      t()
    
    input_df <- cbind(Y, X) %>% as.data.frame()
    colnames(input_df) <- colnames_
    
    model_ <- lm(response ~ ., data = input_df)
    
    print(summary(model_))
    
    model_60days[[type]][[group]] <- model_
    
  }
}


model_30days_data$B$raw_share
model_30days_data$B$raw_abs

summary(model_60days$raw_share$XBB)



