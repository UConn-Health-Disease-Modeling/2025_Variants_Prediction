# this document is to check the smoothing process provided in 'functions_0706'

# set root directory
setwd("/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis")

# load packages 
library(magrittr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(viridis)
library(scales) 
library(dplyr)
library(purrr)

# source the functions (functions_0706.R)
source(glue::glue('/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis/Background/functions_0706.R'))

# select the country
country_select <- 'United Kingdom'
# country_select <- 'South Korea'

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

dataset <- data %>% 
  dplyr::select(country, date, lineage, numerator, unaliased_lineage) %>% 
  dplyr::rename(n = numerator)

# get the classifications 
append_url <- "Data/Aug02_classification_to_append_UnitedKingdom.rds"
classification_to_append <- readRDS(append_url) # 76 classified_label s 

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
  dplyr::filter(classified_label != 'Other') 


# the following result is the sum of monthly share for United Kingdom classified labels  

# # A tibble: 76 × 2
# classified_label               acc_sum
# <chr>                            <dbl>
# 1 Q (B.1.1.7)                     6.11
# 2 JN.1 (BA.2.86.1.1)              5.99
# 3 AY.4 (B.1.617.2.4)              4.72
# 4 B.1.1                           3.94
# 5 XBB.1.5                         3.18
# 6 XBB.1.16                        2.98
# 7 BA.2 (BA.2)                     2.36
# 8 EG.5.1 (XBB.1.9.2.5.1)          2.03
# 9 B.1                             1.77
# 10 BQ.1.1 (BA.5.3.1.1.1.1.1.1)    1.74

# decide the check the smoothing process of the variants: "Q (B.1.1.7)", "JN.1 (BA.2.86.1.1)" and "AY.4 (B.1.617.2.4)"
# Define the variants
# Define the variants
variants <- c("Q (B.1.1.7)", "JN.1 (BA.2.86.1.1)", "AY.4 (B.1.617.2.4)")

# # Filter the classified data and split into datasets for each variant
# variant_datasets <- classified_data %>%
#   dplyr::filter(classified_label %in% variants) %>%
#   split(.$classified_label)
# 
# # Load the smoothed data
# output_url <- "Data/output_UnitedKingdom.rds"
# output <- readRDS(output_url)
# 
# # Filter the output data and split into datasets for each variant
# output_datasets <- output %>%
#   dplyr::filter(classified_label %in% variants) %>%
#   split(.$classified_label)
# 
# # Perform the left join operation for each variant using map
# output_combined <- map(variants, function(variant) {
#   left_join(
#     output_datasets[[variant]],
#     variant_datasets[[variant]],
#     by = "date"
#   )
# })
# 
# # Access the resulting datasets
# output1 <- output_combined[[1]]
# output2 <- output_combined[[2]]
# output3 <- output_combined[[3]]

# so far the numerator is good for classified 'n'
# check the model fitting step by step 




lineages_to_model <- unique(classified_data$classified_label)

# create the output dataframe
output <- data.frame()

# select the lineage 
lin <- variants[1]

# get the lineage data
this_lineage <- classified_data |>
  dplyr::mutate(denominator = sum(n), .by = date) |>
  dplyr::summarise(numerator = sum(n * (classified_label == lin)),
                   .by = c(date, denominator))

# find the min date (set restrictions)
min_date <- this_lineage |>
  dplyr::mutate(
    month = lubridate::floor_date(date, 'month')
  ) |>
  dplyr::mutate(
    total = sum(numerator),
    .by = 'month'
  ) |>
  dplyr::filter(
    numerator >= 1 & total >= 2
  ) |>
  dplyr::filter(
    date == min(date)
  ) |>
  dplyr::pull(
    date
  )

# prepare the data for model 
to_model <- this_lineage |>
  dplyr::filter(
    date >= min_date
  )

# run the model
# output <- growth_rate(
#   numerator = to_model$numerator,
#   denominator = to_model$denominator,
#   date = to_model$date,
#   k_scaling = 10
# ) |>
#   dplyr::mutate(classified_label = lin) |>
#   dplyr::left_join(this_lineage) |>
#   dplyr::bind_rows(output)


numerator = to_model$numerator
denominator = to_model$denominator
date = to_model$date
k_scaling = 10
basis = 'gp'

data <- tidyr::tibble(
  numerator = numerator,
  denominator = denominator,
  date = date
) |>
  dplyr::mutate(
    time_numeric = as.numeric(date - min(date)),
    time_factor = as.factor(time_numeric)
  )

k <- ceiling(max(data$time_numeric) / k_scaling)


formula <- as.formula(
  cbind(numerator, denominator - numerator) ~
    s(time_numeric, k = k, bs = basis)
  + s(time_factor, bs = "re")
)
model <- mgcv::bam(
  formula = formula,
  data = data,
  discrete = T,
  family = "binomial"(
    link = "logit"
  )
)

new_data <- expand.grid(
  date = seq(
    min(date),
    max(date),
    by = 1)
) |>
  dplyr::mutate(
    date = lubridate::ymd(date),
    time_numeric = as.numeric(date - min(date)),
    time_factor = as.factor(time_numeric)
  )

fits <- mgcv::predict.gam(
  model,
  newdata = new_data,
  exclude = "s(time_factor)",
  se.fit = T
)


new_data$fits <- ilogit(fits$fit)
new_data <- new_data %>% 
  dplyr::left_join(
    data, 
    by = "date"
  ) %>% 
  dplyr::mutate(share = ifelse(denominator == 0, NA, numerator / denominator))

# data_sorted <- to_model %>%
#   arrange(date)
ggplot(new_data, aes(x = date)) +
  # Point plot for 'share'
  geom_point(aes(y = share), color = "blue", size = 3) +
  # Line plot for 'fit'
  geom_line(aes(y = fits), color = "red", linewidth = 1) +
  # Add labels and title
  labs(x = "Date", y = "Value", title = "Point Plot for Share and Line Plot for Fit") +
  # Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

################################################################################
############# now understand and directly jump into the data ###################
################################################################################

output_url <- "Data/output_UnitedKingdom.rds"
output <- readRDS(output_url)

lin <- sample(lineages_to_model, 1)

plot_sample <- output %>% 
  dplyr::filter(
    classified_label == lin
  ) %>% 
  dplyr::mutate(
    share = ifelse(denominator == 0, NA, numerator / denominator), 
    fits = ilogit(fit_logit)
  )



ggplot(plot_sample, aes(x = date)) +
  # Point plot for 'share'
  geom_point(aes(y = share), color = "blue", size = 3) +
  # Line plot for 'fit'
  geom_line(aes(y = fits), color = "red", linewidth = 1) +
  # Add labels and title
  labs(x = "Date", y = "Value", title = "Point Plot for Share and Line Plot for Fit") +
  # Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

