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

# set p limit
p_lim <- 0.001

# set n limit
n_lim <- 50

# get the classification to append
# save the result since the function costs a lot of time
classification_to_append <- get_classification(
  lineages = summarised$unaliased_lineage,
  number_sequences = summarised$n,
  p_lim = p_lim, n_lim = n_lim, alias_list
) |>
  dplyr::mutate(
    decimal_lineage = paste0(classified_unasliased, '.')
  )
# 
# saveRDS(classification_to_append, "Data/Jul26_classification_to_append_SouthKorea3.rds")

# load the .rds file 
append_url <- "Data/Jul22_classification_to_append_SouthKorea2.rds"
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

# # Jul22: get 29 'classified_label' for South Korea
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
# saveRDS(output, "Data/output_SouthKorea3.rds")





# load the data
output_url <- "Data/output_SouthKorea3.rds"
output <- readRDS(output_url)

# prefer to remove B, because it accounted 100% share at the start (the initial variant)
# output <- output %>% 
#   dplyr::filter(classified_label != "B")

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
color_palette <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"), brewer.pal(3, "Set1"))

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

dominant_index <- monthly_output %>% 
  dplyr::group_by(classified_label) %>% 
  dplyr::summarise(
    acc_sum = sum(share, na.rm = TRUE), 
    .groups = 'drop'
  ) %>% 
  dplyr::arrange(
    desc(acc_sum)
  ) 

# print(dominant_index, n = 30)
# classified_label           acc_sum
# <chr>                        <dbl>
# 1 B.1.497                      11.4   
# 2 B                            3.59  
# 3 JN.1 (BA.2.86.1.1)           3.58  
# 4 AY.69 (B.1.617.2.69)         3.27  
# 5 BA.5.2 (BA.5.2)              2.90  
# 6 BA.2.3 (BA.2.3)              2.42  
# 7 AY (B.1.617.2)               2.34  
# 8 BF (BA.5.2.1)                2.25  
# 9 HK.3 (XBB.1.9.2.5.1.1.3)     2.21  
# 10 BA.1.1 (BA.1.1)             1.68  
# 11 B.1.619.1                   1.63  
# 12 BA.2 (BA.2)                 1.53  
# 13 JN.1.4 (BA.2.86.1.1.4)      1.50  
# 14 BN.1.2 (BA.2.75.5.1.2)      1.33  
# 15 HK (XBB.1.9.2.5.1.1)        1.32  
# 16 BN.1.3 (BA.2.75.5.1.3)      1.25  
# 17 AY.122.5 (B.1.617.2.122.5)  1.10  
# 18 XBB.1.5                     1.09  
# 19 FL (XBB.1.9.1)              0.950 
# 20 EG.5.1 (XBB.1.9.2.5.1)      0.923 
# 21 XBB.1.16                    0.921 
# 22 EG.1 (XBB.1.9.2.1)          0.809 
# 23 BQ.1 (BA.5.3.1.1.1.1.1)     0.569 
# 24 HG (XBB.2.3.8)              0.493 
# 25 BA.5 (BA.5)                 0.451 
# 26 EG (XBB.1.9.2)              0.314 
# 27 BA.2.75 (BA.2.75)           0.0621
# 28 XBB                         0.0582
# 29 XBB.1.9                     0.0114


# join to the 'output'
output <- output %>% 
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::mutate(
    first_date = min(date, na.rm = TRUE)
  ) %>%
  dplyr::mutate(
    days_since_first = as.integer(date - first_date)
  ) %>%
  dplyr::select(
    -first_date
  ) %>%
  dplyr::ungroup() 

# check the duration of the variants
max_days <- output %>%
  dplyr::group_by(classified_label) %>%
  dplyr::summarise(max_days_since_first = max(days_since_first, na.rm = TRUE), .groups = 'drop')

print(max_days, n = 29)

# # A tibble: 29 × 2
# classified_label           max_days_since_first
# <chr>                                     <int>
# 1 AY (B.1.617.2)                               1111
# 2 AY.122.5 (B.1.617.2.122.5)                   1000
# 3 AY.69 (B.1.617.2.69)                         1085
# 4 B                                            1554
# 5 B.1.497                                      1460
# 6 B.1.619.1                                    1138
# 7 BA.1.1 (BA.1.1)                              880
# 8 BA.2 (BA.2)                                  880
# 9 BA.2.3 (BA.2.3)                              847
# 10 BA.2.75 (BA.2.75)                           667
# 11 BA.5 (BA.5)                                 726
# 12 BA.5.2 (BA.5.2)                             724
# 13 BF (BA.5.2.1)                               725
# 14 BN.1.2 (BA.2.75.5.1.2)                      596
# 15 BN.1.3 (BA.2.75.5.1.3)                      598
# 16 BQ.1 (BA.5.3.1.1.1.1.1)                     619
# 17 EG (XBB.1.9.2)                              451
# 18 EG.1 (XBB.1.9.2.1)                          472
# 19 EG.5.1 (XBB.1.9.2.5.1)                      389
# 20 FL (XBB.1.9.1)                              479
# 21 HG (XBB.2.3.8)                              406
# 22 HK (XBB.1.9.2.5.1.1)                        354
# 23 HK.3 (XBB.1.9.2.5.1.1.3)                    305
# 24 JN.1 (BA.2.86.1.1)                          192
# 25 JN.1.4 (BA.2.86.1.1.4)                      166
# 26 XBB                                         591
# 27 XBB.1.16                                    437
# 28 XBB.1.5                                     519
# 29 XBB.1.9                                     563

# # extract the patterns of first two months (90 days)
# output_90d <- output %>% 
#   dplyr::filter(
#     days_since_first <= 90
#   ) %>% 
#   dplyr::distinct() %>% 
#   dplyr::group_by(
#     classified_label
#   ) %>%
#   slice(-1) %>%
#   ungroup() %>% 
#   dplyr::select(
#     classified_label, date, gr, gr_lower, gr_upper, days_since_first, dominant
#   )
# 
# # visualize (dominant variants (red) v.s. non-dominant variants(blue))
# # Define color scales
# red_colors <- c("yes" = "red")
# blue_colors <- c("no" = "blue")
# 
# # Plot
# ggplot(output_90d, aes(x = days_since_first, y = gr, group = classified_label, color = dominant)) +
#   geom_ribbon(aes(ymin = gr_lower, ymax = gr_upper, fill = dominant), alpha = 0.2, color = NA) +
#   geom_line() +
#   facet_wrap(~ dominant) +
#   scale_color_manual(values = c("yes" = "red", "no" = "blue")) +
#   scale_fill_manual(values = c("yes" = "red", "no" = "blue"), guide = "none") +
#   theme_minimal() +
#   labs(title = "Growth Rate Over Time",
#        x = "Days Since First",
#        y = "Growth Rate",
#        color = "Whether Dominant")


# test_output <- output %>% 
#   dplyr::filter(
#     classified_label == "XBB.1.9"
#   ) %>% 
#   dplyr::select(
#     classified_label, date, prevalence, days_since_first
#   ) 

##########################################
########## correlation analysis ##########
##########################################
# Based on the prevalence of the variant
output_filtered <- output %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    any(prevalence >= 0.05)
  ) # still keeps 26 variants 

# define the initial stage of variants 
output_initial <- output_filtered %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    cummin(prevalence <= 0.05) == 1
  ) %>% 
  dplyr::filter(
    days_since_first != 0
  )

# count the days to reach the 0.05 prevalence threshold 
table(output_initial$classified_label)
# AY.122.5 (B.1.617.2.122.5)       AY.69 (B.1.617.2.69)                  B.1.619.1            BA.1.1 (BA.1.1) 
#                         32                         31                          7                         14 
# BA.2 (BA.2)            BA.2.3 (BA.2.3)                BA.5 (BA.5)            BA.5.2 (BA.5.2) 
#          51                         27                         29                         31 
# BF (BA.5.2.1)     BN.1.2 (BA.2.75.5.1.2)     BN.1.3 (BA.2.75.5.1.3)    BQ.1 (BA.5.3.1.1.1.1.1) 
#            36                         67                         54                         80 
# EG (XBB.1.9.2)         EG.1 (XBB.1.9.2.1)     EG.5.1 (XBB.1.9.2.5.1)             FL (XBB.1.9.1) 
#             66                         68                         68                         58 
# HG (XBB.2.3.8)       HK (XBB.1.9.2.5.1.1)   HK.3 (XBB.1.9.2.5.1.1.3)         JN.1 (BA.2.86.1.1) 
#             41                         41                         45                         33 
# JN.1.4 (BA.2.86.1.1.4)                   XBB.1.16                    XBB.1.5 
#                     19                         42                         69 

# output_initial %>% 
#   dplyr::filter(classified_label == "B.1.619.1")


# hypothesis1: the time to the threshold (0.05)
time2threshold <- output_initial %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::mutate(
    freq = dplyr::n()
  ) %>% 
  dplyr::select(classified_label, freq) %>% 
  dplyr::distinct() %>% 
  dplyr::left_join(dominant_index, by = "classified_label")
  

cor_test_result <- cor.test(time2threshold$freq, time2threshold$acc_sum)
print(cor_test_result)

# Pearson's product-moment correlation
# 
# data:  time2threshold$freq and time2threshold$acc_sum
# t = -2.468, df = 21, p-value = 0.02226
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.7414482 -0.0770127
# sample estimates:
#        cor 
# -0.4741627 


# hypothesis2: average growth rate at initial stage 
avg_gr <- output_initial %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::mutate(
    avg_gr = mean(gr)
  ) %>% 
  dplyr::select(classified_label, avg_gr) %>% 
  dplyr::distinct() %>% 
  dplyr::left_join(dominant_index, by = "classified_label")

cor_test_result <- cor.test(avg_gr$avg_gr, avg_gr$acc_sum)
print(cor_test_result)

# Pearson's product-moment correlation
# 
# data:  avg_gr$avg_gr and avg_gr$acc_sum
# t = 1.2608, df = 21, p-value = 0.2212
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1649610  0.6107021
# sample estimates:
#       cor 
# 0.2652794 
