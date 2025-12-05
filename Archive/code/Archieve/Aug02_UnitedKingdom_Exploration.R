# source the functions (functions_0706.R)
source(glue::glue('/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis/Background/functions_0706.R'))

# load packages 
library(magrittr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(viridis)
library(scales) 

# select the country
country_select <- 'United Kingdom'

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

# summarised <- dataset %>% 
#   # dplyr::filter(date >= (max(date) - date_range)) |>
#   dplyr::summarise(n = sum(n), .by = unaliased_lineage)

# # set p limit
# p_lim <- 0.001
# 
# # set n limit
# n_lim <- 50

# get the classification to append
# save the result since the function costs a lot of time
# classification_to_append <- get_classification(
#   lineages = summarised$unaliased_lineage,
#   number_sequences = summarised$n,
#   p_lim = p_lim, n_lim = n_lim, alias_list
# ) |>
#   dplyr::mutate(
#     decimal_lineage = paste0(classified_unasliased, '.')
#   )

# saveRDS(classification_to_append, "Data/Aug02_classification_to_append_UnitedKingdom.rds")

# # load the .rds file
# append_url <- "Data/Aug02_classification_to_append_UnitedKingdom.rds"
# classification_to_append <- readRDS(append_url)
# 
# # append the classification
# classified_data <- dataset |>
#   dplyr::mutate(
#     decimal_lineage = paste0(unaliased_lineage, '.')
#   ) |>
#   fuzzyjoin::fuzzy_left_join(
#     classification_to_append,
#     by = 'decimal_lineage',
#     match_fun = stringr::str_starts
#   ) |>
#   dplyr::filter(
#     is.na(classified_unasliased) | stringr::str_length(classified_unasliased) == max(stringr::str_length(classified_unasliased)),
#     .by = c(lineage, date)
#   ) |>
#   dplyr::mutate(
#     length_class = stringr::str_count(classified_unasliased, '\\.'),
#     length_lineage = stringr::str_count(unaliased_lineage, '\\.'),
#     classified_label = dplyr::if_else((length_lineage - length_class) <= 1,
#                                       classified_label,
#                                       NA),
#     classified_label = dplyr::if_else(is.na(classified_label),
#                                       "Other",
#                                       classified_label)
#   ) |>
#   dplyr::filter(classified_label != 'Other')
# 
# # Aug02: get 29 'classified_label' for United Kingdom
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

# # save the output for United Kingdom
# saveRDS(output, "Data/output_UnitedKingdom.rds")







# load the data
output_url <- "Data/output_UnitedKingdom.rds"
output <- readRDS(output_url)

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
colors_brewer <- c(
  brewer.pal(12, "Set3"),  # 12 colors
  brewer.pal(8, "Dark2"),  # 8 colors
  brewer.pal(8, "Set2"),   # 8 colors
  brewer.pal(9, "Set1"),   # 9 colors (max available in Set1)
  brewer.pal(8, "Pastel1"),# 8 colors
  brewer.pal(8, "Pastel2"),# 8 colors
  brewer.pal(11, "Spectral") # 11 colors
)

# Check the current number of colors
length(colors_brewer)  # Total = 64 colors

# Generate additional colors to reach 76 using colorRampPalette
additional_colors <- colorRampPalette(c("blue", "red", "green", "orange", "purple", "brown"))(12)

# Combine to form the final color palette of 76 colors
color_palette <- c(colors_brewer, additional_colors)

ggplot(monthly_output, aes(x = month, y = share, fill = classified_label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "Monthly Share of Classified Labels",
    x = "Month",
    y = "Share",
    fill = "Classified Label"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"  # This line removes the legend
  )

dominant_index <- monthly_output %>% 
  dplyr::group_by(classified_label) %>% 
  dplyr::summarise(
    acc_sum = sum(share, na.rm = TRUE), 
    .groups = 'drop'
  ) %>% 
  dplyr::arrange(
    desc(acc_sum)
  ) 

print(dominant_index)

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


