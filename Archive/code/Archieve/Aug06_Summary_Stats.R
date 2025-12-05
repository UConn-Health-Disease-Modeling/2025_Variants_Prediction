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
# threshold 0.05

Korea_prev <- Korea_data %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    max(sharing) >= 0.05
  ) %>%
  dplyr::mutate(
    first_sharing_date = min(date[sharing >= 0.05])
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    date <= first_sharing_date
  ) %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::summarise(
    appearances = n(),
    lowest_sharing = min(sharing),
    highest_sharing = max(sharing)
  ) %>% 
  dplyr::left_join(
    Korea_dominance, 
    by = "classified_label"
  )

UK_prev <- UK_data %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    max(sharing) >= 0.05
  ) %>%
  dplyr::mutate(
    first_sharing_date = min(date[sharing >= 0.05])
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    date <= first_sharing_date
  ) %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::summarise(
    appearances = n(),
    lowest_sharing = min(sharing),
    highest_sharing = max(sharing)
  ) %>% 
  dplyr::left_join(
    UK_dominance, 
    by = "classified_label"
  )

Korea_model <- lm(dominance ~ appearances, data = Korea_prev)
UK_model <- lm(dominance ~ appearances, data = UK_prev)
summary(Korea_model)
summary(UK_model)
################################################################################
# threshold 0.1

Korea_prev <- Korea_data %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    max(sharing) >= 0.1
  ) %>%
  dplyr::mutate(
    first_sharing_date = min(date[sharing >= 0.1])
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    date <= first_sharing_date
  ) %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::summarise(
    appearances = n(),
    lowest_sharing = min(sharing),
    highest_sharing = max(sharing)
  ) %>% 
  dplyr::left_join(
    Korea_dominance, 
    by = "classified_label"
  )

UK_prev <- UK_data %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    max(sharing) >= 0.1
  ) %>%
  dplyr::mutate(
    first_sharing_date = min(date[sharing >= 0.1])
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    date <= first_sharing_date
  ) %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::summarise(
    appearances = n(),
    lowest_sharing = min(sharing),
    highest_sharing = max(sharing)
  ) %>% 
  dplyr::left_join(
    UK_dominance, 
    by = "classified_label"
  )

Korea_model <- lm(dominance ~ appearances, data = Korea_prev)
UK_model <- lm(dominance ~ appearances, data = UK_prev)
summary(Korea_model)
summary(UK_model)