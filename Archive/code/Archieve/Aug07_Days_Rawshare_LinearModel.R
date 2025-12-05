# source the functions (functions_0706.R)
source(glue::glue('/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis/Background/functions_0706.R'))
source('Code/functions_franky_0725.R')

# load packages 
library(magrittr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(viridis)
library(scales) 
library(fda)

################################################################################
## Classified the data (for South Koera and United Kingdom together)

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


# # classified the variants
# summarised <- dataset %>%
#   dplyr::group_by(country, unaliased_lineage) %>% 
#   dplyr::summarise(n = sum(n))
# 
# # set p limit
# p_lim <- 0.001
# 
# # set n limit
# n_lim <- 50
# 
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
# saveRDS(classification_to_append, "Data/Aug06_classification_to_append_Combined.rds")

# load the .rds file
append_url <- "Data/Aug06_classification_to_append_Combined.rds"
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
  dplyr::select(country, date, lineage, n, classified_label)

# check the number of the classified labels
# length(unique(classified_data$lineage)) # 1608
# length(unique(classified_data$classified_label)) # 88

classified_data <- classified_data %>% 
  dplyr::select(
    country, date, n, classified_label
  ) %>% 
  dplyr::group_by(
    country, classified_label, date, 
  ) %>% 
  dplyr::summarise(
    n_total = sum (n)
  )


head(classified_data)
colnames(classified_data)

# complete the dates to be continuous
complete_dates <- function(data) {
  data %>%
    dplyr::group_by(
      country, 
      classified_label
    ) %>%
    tidyr::complete(date = seq.Date(min(date), max(date), by = "day")
    ) %>%
    tidyr::replace_na(
      list(n_total = 0)
    ) %>%
    dplyr::ungroup()
}

classified_data <- complete_dates(classified_data) %>% 
  dplyr::rename(
    numerator = n_total
  ) %>% 
  dplyr::group_by(
    country, date
  ) %>% 
  dplyr::mutate(
    denominator = sum(numerator, na.rm = TRUE)
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    share = numerator/denominator
  )

# separate data
Korea_data <- classified_data %>%
  dplyr::filter(
    country == "South Korea"
  ) %>%
  dplyr::arrange(
    classified_label, date
  ) %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    dplyr::n() > 30
  ) %>%
  dplyr::mutate(
    sharing = dplyr::coalesce(share, 0),
    gap = as.numeric(date - dplyr::first(date))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-share)

Korea_dominance <- Korea_data %>% 
  dplyr::group_by(
    classified_label
  ) %>% 
  dplyr::summarise(
    dominance = sum(sharing, na.rm = TRUE),
  ) %>% 
  dplyr::arrange(
    dominance
  )

UK_data <- classified_data %>%
  dplyr::filter(
    country == "United Kingdom"
  ) %>%
  dplyr::arrange(
    classified_label, date
  ) %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    dplyr::n() > 30
  ) %>%
  dplyr::mutate(
    sharing = dplyr::coalesce(share, 0),
    gap = as.numeric(date - dplyr::first(date))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-share)

UK_dominance <- UK_data %>% 
  dplyr::group_by(
    classified_label
  ) %>% 
  dplyr::summarise(
    dominance = sum(sharing, na.rm = TRUE),
  ) %>% 
  dplyr::arrange(
    dominance
  )

# length(unique(Korea_data$classified_label)) # 80
# length(unique(UK_data$classified_label)) # 86

# check if the share is good (Bingo!)
# sample_date <- sample(unique(Korea_data$date), 5)
# sample_index <- sample(1:length(sample_date), 1)
# sum((Korea_data %>% dplyr::filter(date == sample_date[sample_index]))$share)

rm(alias_list, classified_data, data, append_url, url)





################################################################################
################################################################################
# 30 DAYS

# For Korea data: 
korea_list <- fda_matrix_prepare(days_input = 30, Korea_data, Korea_dominance)
Korea_names <- korea_list$names
Korea_input <- korea_list$input
Korea_response <- korea_list$response[1:length(Korea_names), ]

# For UK data:
uk_list <- fda_matrix_prepare(days_input = 30, UK_data, UK_dominance)
UK_names <- uk_list$names
UK_input <- uk_list$input
UK_response <- uk_list$response[1:length(UK_names), ]

# grab the gr (growth rates)
UK_output <- BSpline_gr(UK_input)
UK_fitted <- UK_output$fitted
UK_gr <- UK_output$gr
UK_raw <- UK_output$raw

Korea_output <- BSpline_gr(Korea_input)
Korea_fitted <- Korea_output$fitted
Korea_gr <- Korea_output$gr
Korea_raw <- Korea_output$raw

Korea_compress_share <- apply(Korea_raw, 2, paa_transform, num_segments = 6)
UK_compress_share    <- apply(UK_raw,    2, paa_transform, num_segments = 6)

## Do the Linear regression model 
Korea_df <- as.data.frame(t(Korea_compress_share)) %>% 
  dplyr::rename(
    '1-5 share' = V1, 
    '6-10 share' = V2, 
    '11-15 share' = V3, 
    '16-20 share' = V4, 
    '21-25 share' = V5, 
    '26-30 share' = V6
  )
Korea_df$dominance <- Korea_response

UK_df <- as.data.frame(t(UK_compress_share)) %>% 
  dplyr::rename(
    '1-5 share' = V1, 
    '6-10 share' = V2, 
    '11-15 share' = V3, 
    '16-20 share' = V4, 
    '21-25 share' = V5, 
    '26-30 share' = V6
  )
UK_df$dominance <- UK_response

Korea_model <- lm(dominance ~ ., data = Korea_df)
summary(Korea_model)

UK_model <- lm(dominance ~ ., data = UK_df)
summary(UK_model)

Korea_corr <- correlation_tests(Korea_df, "dominance")
UK_corr <- correlation_tests(UK_df, "dominance")





################################################################################
################################################################################
# 45 DAYS

# For Korea data: 
korea_list <- fda_matrix_prepare(days_input = 45, Korea_data, Korea_dominance)
Korea_names <- korea_list$names
Korea_input <- korea_list$input
Korea_response <- korea_list$response[1:length(Korea_names), ]

# For UK data:
uk_list <- fda_matrix_prepare(days_input = 45, UK_data, UK_dominance)
UK_names <- uk_list$names
UK_input <- uk_list$input
UK_response <- uk_list$response[1:length(UK_names), ]

# grab the gr (growth rates)
UK_output <- BSpline_gr(UK_input)
UK_fitted <- UK_output$fitted
UK_gr <- UK_output$gr
UK_raw <- UK_output$raw

Korea_output <- BSpline_gr(Korea_input)
Korea_fitted <- Korea_output$fitted
Korea_gr <- Korea_output$gr
Korea_raw <- Korea_output$raw


Korea_compress_share <- apply(Korea_raw, 2, paa_transform, num_segments = 9)
UK_compress_share    <- apply(UK_raw,    2, paa_transform, num_segments = 9)

## Do the Linear regression model 
Korea_df <- as.data.frame(t(Korea_compress_share)) %>% 
  dplyr::rename(
    '1-5 share' = V1, 
    '6-10 share' = V2, 
    '11-15 share' = V3, 
    '16-20 share' = V4, 
    '21-25 share' = V5, 
    '26-30 share' = V6, 
    '31-35 share' = V7, 
    '36-40 share' = V8, 
    '41-45 share' = V9
  )
Korea_df$dominance <- Korea_response

UK_df <- as.data.frame(t(UK_compress_share)) %>% 
  dplyr::rename(
    '1-5 share' = V1, 
    '6-10 share' = V2, 
    '11-15 share' = V3, 
    '16-20 share' = V4, 
    '21-25 share' = V5, 
    '26-30 share' = V6, 
    '31-35 share' = V7, 
    '36-40 share' = V8, 
    '41-45 share' = V9
  )
UK_df$dominance <- UK_response

Korea_model <- lm(dominance ~ ., data = Korea_df)
summary(Korea_model)

UK_model <- lm(dominance ~ ., data = UK_df)
summary(UK_model)

Korea_corr <- correlation_tests(Korea_df, "dominance")
UK_corr <- correlation_tests(UK_df, "dominance")





################################################################################
################################################################################
# 60 DAYS

# For Korea data: 
korea_list <- fda_matrix_prepare(days_input = 60, Korea_data, Korea_dominance)
Korea_names <- korea_list$names
Korea_input <- korea_list$input
Korea_response <- korea_list$response[1:length(Korea_names), ]

# For UK data:
uk_list <- fda_matrix_prepare(days_input = 60, UK_data, UK_dominance)
UK_names <- uk_list$names
UK_input <- uk_list$input
UK_response <- uk_list$response[1:length(UK_names), ]

# grab the gr (growth rates)
UK_output <- BSpline_gr(UK_input)
UK_fitted <- UK_output$fitted
UK_gr <- UK_output$gr
UK_raw <- UK_output$raw

Korea_output <- BSpline_gr(Korea_input)
Korea_fitted <- Korea_output$fitted
Korea_gr <- Korea_output$gr
Korea_raw <- Korea_output$raw


Korea_compress_share <- apply(Korea_raw, 2, paa_transform, num_segments = 12)
UK_compress_share    <- apply(UK_raw,    2, paa_transform, num_segments = 12)

## Do the Linear regression model 
Korea_df <- as.data.frame(t(Korea_compress_share)) %>% 
  dplyr::rename(
    '1-5 share' = V1, 
    '6-10 share' = V2, 
    '11-15 share' = V3, 
    '16-20 share' = V4, 
    '21-25 share' = V5, 
    '26-30 share' = V6, 
    '31-35 share' = V7, 
    '36-40 share' = V8, 
    '41-45 share' = V9, 
    '46-50 share' = V10, 
    '51-55 share' = V11, 
    '56-60 share' = V12
  )
Korea_df$dominance <- Korea_response

UK_df <- as.data.frame(t(UK_compress_share)) %>% 
  dplyr::rename(
    '1-5 share' = V1, 
    '6-10 share' = V2, 
    '11-15 share' = V3, 
    '16-20 share' = V4, 
    '21-25 share' = V5, 
    '26-30 share' = V6, 
    '31-35 share' = V7, 
    '36-40 share' = V8, 
    '41-45 share' = V9, 
    '46-50 share' = V10, 
    '51-55 share' = V11, 
    '56-60 share' = V12
  )
UK_df$dominance <- UK_response

Korea_model <- lm(dominance ~ ., data = Korea_df)
summary(Korea_model)

UK_model <- lm(dominance ~ ., data = UK_df)
summary(UK_model)

Korea_corr <- correlation_tests(Korea_df, "dominance")
UK_corr <- correlation_tests(UK_df, "dominance")
