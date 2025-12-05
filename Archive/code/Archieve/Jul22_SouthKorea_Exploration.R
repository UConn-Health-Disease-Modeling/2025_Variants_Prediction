# source the functions (functions_0706.R)
source(glue::glue('/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis/Background/functions_0706.R'))

# set root directory
setwd("/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis")

# load packages 
library(magrittr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(viridis)

# select the country
country_select <- 'South Korea'

# get the alias
alias_list <- get_alias()

# set the time range (as large as possible)
months_from_most_recent <- 300 * 7

# load the data and process 
url <- 'Data/summary_gisaid_new.csv'
data <- readr::read_csv(url) |>
  dplyr::filter(
    country == country_select,
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

# # check the data 
# unique(data$country)

dataset <- data %>% 
  dplyr::select(country, date, lineage, numerator, unaliased_lineage) %>% 
  dplyr::rename(n = numerator)

# no time range needed to be set here 
# date_range <- 300 * 7

# summaries the dataset
summarised <- dataset |>
  # dplyr::filter(date >= (max(date) - date_range)) |>
  dplyr::summarise(n = sum(n), .by = unaliased_lineage)

# # set p limit
# p_lim <- 0.03
# 
# # set n limit
# n_lim <- 1000
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
# saveRDS(classification_to_append, "Data/Jul22_classification_to_append_SouthKorea.rds")

append_url <- "Data/Jul22_classification_to_append_SouthKorea.rds"
classification_to_append <- readRDS(append_url)
# unique(classification_to_append$classified_label)

# process the data with 'classification_to_append' (Jul22: 32,018 obs to 10,860)
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
  dplyr::filter(classified_label != 'Other')

# # Jul22: get 14 'classified_label' for South Korea
# lineages_to_model <- unique(classified_data$classified_label)
# 
# # create the output dataframe
# output <- data.frame()
# 
# # this part of code is already Given
# for(lin in lineages_to_model){
# 
#   this_lineage <- classified_data |>
#     dplyr::mutate(denominator = sum(n), .by = date) |>
#     dplyr::summarise(numerator = sum(n * (classified_label == lin)),
#                      .by = c(date, denominator))
# 
#   min_date <- this_lineage |>
#     dplyr::mutate(
#       month = lubridate::floor_date(date, 'month')
#     ) |>
#     dplyr::mutate(
#       total = sum(numerator),
#       .by = 'month'
#     ) |>
#     dplyr::filter(
#       numerator >= 1 & total >= 2
#     ) |>
#     dplyr::filter(
#       date == min(date)
#     ) |>
#     dplyr::pull(
#       date
#     )
# 
#   to_model <- this_lineage |>
#     dplyr::filter(
#       date >= min_date
#     )
# 
#   output <- growth_rate(
#     numerator = to_model$numerator,
#     denominator = to_model$denominator,
#     date = to_model$date,
#     k_scaling = 10
#   ) |>
#     dplyr::mutate(classified_label = lin) |>
#     dplyr::left_join(this_lineage) |>
#     dplyr::bind_rows(output)
# 
# 
# }
# 
# # save the output for South Korea
# saveRDS(output, "Data/output_SouthKorea.rds")





# load the data
output_url <- "Data/output_SouthKorea.rds"
output <- readRDS(output_url)

# prefer to remove B, because it accounted 100% share at the start (the initial variant)
output <- output %>%
dplyr::filter(classified_label != "B")


## calculate the monthly share
monthly_output <-output %>%
  dplyr::mutate(
    month = floor_date(date, "month")
  ) %>%
  dplyr::group_by(
    classified_label, month
  ) %>%
  dplyr::  summarise(
    monthly_sum = sum(numerator, na.rm = TRUE), 
    .groups = 'drop'
  )

# Calculate total monthly sum
monthly_output_sum <- monthly_output %>%
  dplyr::group_by(month) %>%
  dplyr::summarise(total_sum = sum(monthly_sum, na.rm = TRUE), .groups = 'drop')

# Merge the monthly sums with the total monthly sums
monthly_output <- monthly_output %>%
  dplyr::left_join(
    monthly_output_sum, by = "month"
  ) %>%
  dplyr::mutate(
    share = monthly_sum / total_sum
  ) %>%
  dplyr::select(
    classified_label, month, monthly_sum, share
  )


## Create the plot
# Define a color palette with more distinguishable colors
color_palette <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"))

ggplot(monthly_output, aes(x = month, y = share, fill = classified_label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = color_palette) +
  labs(title = "Monthly Share of Classified Labels",
       x = "Month",
       y = "Share",
       fill = "Classified Label") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# check the data again (all good)
# sum((monthly_output %>% dplyr::filter(month == "2022-03-01"))$share)


## calculate the unit growth rate 
# rank the sum of share to define whether dominant or not 
# 
# classified_label         acc_sum
# <chr>                      <dbl>
# 1 JN.1 (BA.2.86.1.1)         5.18 
# 2 AY (B.1.617.2)             4.79 
# 3 AY.69 (B.1.617.2.69)       4.29 
# 4 HK.3 (XBB.1.9.2.5.1.1.3)   3.92 
# 5 XBB.1                      3.66 
# 6 BA.5.2 (BA.5.2)            3.20 
# --------------------------------------
# 7 BN.1 (BA.2.75.5.1)         2.71 
# 8 BF (BA.5.2.1)              2.46 
# 9 BA.2.3 (BA.2.3)            2.45 
# 10 BA.1.1 (BA.1.1)            1.90 
# 11 BA.2 (BA.2)                1.57 
# 12 XBB.1.9                    1.40 
# 13 BA.5 (BA.5)                0.466
# 
dominant_index <- monthly_output %>% 
  dplyr::group_by(classified_label) %>% 
  dplyr::summarise(
    acc_sum = sum(share, na.rm = TRUE), 
    .groups = 'drop'
  ) %>% 
  dplyr::arrange(
    desc(acc_sum)
  ) %>% 
  dplyr::mutate(
    dominant = ifelse(acc_sum > 3, "yes", "no")
  ) %>% 
  dplyr::select(
    classified_label, dominant
  )

# join to the 'output'
output <- output %>% 
  dplyr::left_join(
    dominant_index, 
    by = "classified_label"
  ) %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::mutate(
    first_date = min(date, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    days_since_first = as.integer(date - first_date)
  ) %>%
  dplyr::select(-first_date
  ) %>%
  dplyr::ungroup()
  
# # check the duration of the variants 
# max_days <- output %>%
#   group_by(classified_label) %>%
#   summarise(max_days_since_first = max(days_since_first, na.rm = TRUE), .groups = 'drop')
# 
# print(max_days)
# classified_label         max_days_since_first
# <chr>                                   <int>
#   1 AY (B.1.617.2)                        1111
# 2 AY.69 (B.1.617.2.69)                    1085
# 3 BA.1.1 (BA.1.1)                         880
# 4 BA.2 (BA.2)                             880
# 5 BA.2.3 (BA.2.3)                         847
# 6 BA.5 (BA.5)                             726
# 7 BA.5.2 (BA.5.2)                         724
# 8 BF (BA.5.2.1)                           725
# 9 BN.1 (BA.2.75.5.1)                      598
# 10 HK.3 (XBB.1.9.2.5.1.1.3)               305
# 11 JN.1 (BA.2.86.1.1)                     192
# 12 XBB.1                                  591
# 13 XBB.1.9                                563

# extract the patterns of first two months (90 days)
output_90d <- output %>% 
  dplyr::filter(
    days_since_first <= 90
  ) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(
    classified_label
  ) %>%
  slice(-1) %>%
  ungroup() %>% 
  dplyr::select(
    classified_label, date, gr, gr_lower, gr_upper, days_since_first, dominant
  )

# check <- output_60d %>% dplyr::filter(classified_label == "XBB.1.9")

# visualize (dominant variants (red) v.s. non-dominant variants(blue))
# Define color scales
red_colors <- c("yes" = "red")
blue_colors <- c("no" = "blue")

# Plot
ggplot(output_90d, aes(x = days_since_first, y = gr, group = classified_label, color = dominant)) +
  geom_ribbon(aes(ymin = gr_lower, ymax = gr_upper, fill = dominant), alpha = 0.2, color = NA) +
  geom_line() +
  facet_wrap(~ dominant) +
  scale_color_manual(values = c("yes" = "red", "no" = "blue")) +
  scale_fill_manual(values = c("yes" = "red", "no" = "blue"), guide = "none") +
  theme_minimal() +
  labs(title = "Growth Rate Over Time",
       x = "Days Since First",
       y = "Growth Rate",
       color = "Whether Dominant")

# find the outlier for 'yes' dominant variants
# output_60d %>% 
#   dplyr::filter(
#     dominant == "yes", 
#     gr_upper < 0
#   )
# 
# AY (B.1.617.2)














