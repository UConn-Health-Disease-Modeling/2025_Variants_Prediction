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
    desc(dominance)
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
    desc(dominance)
  )

# length(unique(Korea_data$classified_label)) # 80
# length(unique(UK_data$classified_label)) # 86

# check if the share is good (Bingo!)
# sample_date <- sample(unique(Korea_data$date), 5)
# sample_index <- sample(1:length(sample_date), 1)
# sum((Korea_data %>% dplyr::filter(date == sample_date[sample_index]))$share)

rm(alias_list, classified_data, data, append_url, url)





################################################################################
# do the boxplot to check the distribution of variants
process_monthly_data <- function(data) {
  data %>%
    dplyr::mutate(
      month = floor_date(date, "month")
    ) %>%
    dplyr::group_by(classified_label, month) %>%
    dplyr::summarise(
      monthly_sum = sum(numerator, na.rm = TRUE), 
      .groups = 'drop'
    ) %>%
    dplyr::group_by(month) %>%
    dplyr::mutate(
      total_sum = sum(monthly_sum, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      share = monthly_sum / total_sum
    ) %>%
    dplyr::select(
      classified_label, month, monthly_sum, share
    )
}

# Apply the function to both Korea and UK datasets
Korea_monthly <- process_monthly_data(Korea_data)
UK_monthly <- process_monthly_data(UK_data)

Korea_legend <- Korea_dominance$classified_label[1:10]
UK_legend <- UK_dominance$classified_label[1:10]


color_palette_72 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#00CED1", 
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9",
  "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928",
  "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", "#FFFF99", "#FF69B4", "#FF1493", "#FF6347",
  "#FFD700", "#ADFF2F", "#00FF7F", "#00FA9A", "#7FFFD4", "#4682B4", "#6495ED", "#7B68EE", "#9370DB",
  "#8A2BE2", "#BA55D3", "#D8BFD8", "#DDA0DD", "#EE82EE", "#FFC0CB", "#FFB6C1", "#FFA07A", "#FF4500",
  "#DA70D6", "#EE82EE", "#FFD700", "#7FFF00", "#00FF00", "#32CD32", "#FF00FF", "#800080", "#8B0000"
)

ggplot(Korea_monthly, aes(x = month, y = share, fill = classified_label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = color_palette_72) +
  labs(
    title = "Monthly Share (South Korea)",
    x = "Month",
    y = "Share",
    fill = "Classified Label"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.key.size = unit(0.15, "cm"),
    legend.text = element_text(size = 4), 
    legend.title = element_text(size = 4),
    plot.title = element_text(hjust = 0.5, face = "bold") 
  ) +
  guides(fill = guide_legend(nrow = 3))

color_palette_86 <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#00CED1", 
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9",
  "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928",
  "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", "#FFFF99", "#FF69B4", "#FF1493", "#FF6347",
  "#FFD700", "#ADFF2F", "#00FF7F", "#00FA9A", "#7FFFD4", "#4682B4", "#6495ED", "#7B68EE", "#9370DB",
  "#8A2BE2", "#BA55D3", "#D8BFD8", "#DDA0DD", "#EE82EE", "#FFC0CB", "#FFB6C1", "#FFA07A", "#FF4500",
  "#DA70D6", "#EE82EE", "#FFD700", "#7FFF00", "#00FF00", "#32CD32", "#FF00FF", "#800080", "#8B0000",
  "#708090", "#FF8C00", "#CD5C5C", "#4B0082", "#2E8B57", "#00CED1", "#FF4500", "#FA8072", "#20B2AA",
  "#00BFFF", "#FF6347", "#D2691E", "#8B4513", "#4682B4", "#9ACD32"
)

ggplot(UK_monthly, aes(x = month, y = share, fill = classified_label)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = color_palette_86) +
  labs(
    title = "Monthly Share (United Kingdom)",
    x = "Month",
    y = "Share",
    fill = "Classified Label"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.key.size = unit(0.15, "cm"),
    legend.text = element_text(size = 4), 
    legend.title = element_text(size = 4),
    plot.title = element_text(hjust = 0.5, face = "bold") 
  ) +
  guides(fill = guide_legend(nrow = 3))




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








# visualize
plot_share <- function(response, data, title) {
  Korea_dom_index <- response > 45
  
  Korea_fitted_dom <- as.data.frame(data[, Korea_dom_index])
  Korea_fitted_nondom <- as.data.frame(data[, !Korea_dom_index])
  
  Korea_fitted_dom$time <- 1:nrow(Korea_fitted_dom)
  Korea_fitted_nondom$time <- 1:nrow(Korea_fitted_nondom)
  
  Korea_fitted_dom_melt <- reshape2::melt(Korea_fitted_dom, id.vars = "time")
  Korea_fitted_nondom_melt <- reshape2::melt(Korea_fitted_nondom, id.vars = "time")
  
  Korea_fitted_dom_melt$dominance <- "yes"
  Korea_fitted_nondom_melt$dominance <- "no"
  
  plot_data <- rbind(Korea_fitted_dom_melt, Korea_fitted_nondom_melt)
  
  ggplot() +
    geom_line(data = plot_data, aes(x = time, y = value, group = variable, color = dominance), alpha = 0.7) +
    facet_wrap(~dominance) +
    labs(
      title = title,
      x = "Time",
      y = "Value"
    ) +
    theme_bw()
}

plot_share(Korea_response, Korea_fitted, "Fitted Share (South Korea)")
plot_share(UK_response, UK_fitted, "Fitted Share (United Kingdom)")

plot_share(Korea_response, Korea_gr, "Fitted Growth Rate (South Korea)")
plot_share(UK_response, UK_gr, "Fitted Growth Rate (United Kingdom)")






Korea_compress_gr <- apply(Korea_gr, 2, paa_transform, num_segments = 6)
UK_compress_gr    <- apply(UK_gr,    2, paa_transform, num_segments = 6)
# index <- sample(1:50, 1)
# selected_row <- Korea_gr[ , index]
# compressed_row <- Korea_compress_gr[ , index]
# 
# # prepare the data
# original_df <- data.frame(
#   Time = 1:length(selected_row),
#   Value = selected_row,
#   Type = "Original"
# )
# 
# 
# compressed_df <- data.frame(
#   Time = seq(1, length(selected_row), length.out = length(compressed_row)),
#   Value = compressed_row,
#   Type = "Compressed"
# )
# 
# # draw the plot
# ggplot() +
#   geom_line(data = original_df, aes(x = Time, y = Value, color = Type), size = 1) +
#   geom_line(data = compressed_df, aes(x = Time, y = Value, color = Type), size = 1, linetype = "dashed") +
#   labs(title = "Comparison of Original and Compressed Time Series",
#        x = "Time",
#        y = "Value") +
#   scale_color_manual(values = c("Original" = "blue", "Compressed" = "red")) +
#   theme_minimal()

## Do the Linear regression model 
Korea_df <- as.data.frame(t(Korea_compress_gr)) %>% 
  dplyr::rename(
    '1-5 gr' = V1, 
    '6-10 gr' = V2, 
    '11-15 gr' = V3, 
    '16-20 gr' = V4, 
    '21-25 gr' = V5, 
    '26-30 gr' = V6
  )
Korea_df$dominance <- Korea_response

UK_df <- as.data.frame(t(UK_compress_gr)) %>% 
  dplyr::rename(
    '1-5 gr' = V1, 
    '6-10 gr' = V2, 
    '11-15 gr' = V3, 
    '16-20 gr' = V4, 
    '21-25 gr' = V5, 
    '26-30 gr' = V6
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


Korea_compress_gr <- apply(Korea_gr, 2, paa_transform, num_segments = 9)
UK_compress_gr    <- apply(UK_gr,    2, paa_transform, num_segments = 9)

## Do the Linear regression model 
Korea_df <- as.data.frame(t(Korea_compress_gr)) %>% 
  dplyr::rename(
    '1-5 gr' = V1, 
    '6-10 gr' = V2, 
    '11-15 gr' = V3, 
    '16-20 gr' = V4, 
    '21-25 gr' = V5, 
    '26-30 gr' = V6, 
    '31-35 gr' = V7, 
    '36-40 gr' = V8, 
    '41-45 gr' = V9
  )
Korea_df$dominance <- Korea_response

UK_df <- as.data.frame(t(UK_compress_gr)) %>% 
  dplyr::rename(
    '1-5 gr' = V1, 
    '6-10 gr' = V2, 
    '11-15 gr' = V3, 
    '16-20 gr' = V4, 
    '21-25 gr' = V5, 
    '26-30 gr' = V6, 
    '31-35 gr' = V7, 
    '36-40 gr' = V8, 
    '41-45 gr' = V9
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


Korea_compress_gr <- apply(Korea_gr, 2, paa_transform, num_segments = 12)
UK_compress_gr    <- apply(UK_gr,    2, paa_transform, num_segments = 12)

## Do the Linear regression model 
Korea_df <- as.data.frame(t(Korea_compress_gr)) %>% 
  dplyr::rename(
    '1-5 gr' = V1, 
    '6-10 gr' = V2, 
    '11-15 gr' = V3, 
    '16-20 gr' = V4, 
    '21-25 gr' = V5, 
    '26-30 gr' = V6, 
    '31-35 gr' = V7, 
    '36-40 gr' = V8, 
    '41-45 gr' = V9, 
    '46-50 gr' = V10, 
    '51-55 gr' = V11, 
    '56-60 gr' = V12
  )
Korea_df$dominance <- Korea_response

UK_df <- as.data.frame(t(UK_compress_gr)) %>% 
  dplyr::rename(
    '1-5 gr' = V1, 
    '6-10 gr' = V2, 
    '11-15 gr' = V3, 
    '16-20 gr' = V4, 
    '21-25 gr' = V5, 
    '26-30 gr' = V6, 
    '31-35 gr' = V7, 
    '36-40 gr' = V8, 
    '41-45 gr' = V9, 
    '46-50 gr' = V10, 
    '51-55 gr' = V11, 
    '56-60 gr' = V12
  )
UK_df$dominance <- UK_response

Korea_model <- lm(dominance ~ ., data = Korea_df)
summary(Korea_model)

UK_model <- lm(dominance ~ ., data = UK_df)
summary(UK_model)

Korea_corr <- correlation_tests(Korea_df, "dominance")
UK_corr <- correlation_tests(UK_df, "dominance")


# because of the correlations including more days does not necessarily improve the prediction efficiency
