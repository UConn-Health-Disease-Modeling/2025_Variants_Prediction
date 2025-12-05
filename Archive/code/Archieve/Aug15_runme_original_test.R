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

################################################################################
## apply the 'extract_growth_rates' function
type = "absolute-growth" # "relative-growth" or "absolute-growth"

denominator = 10

dow = 'none'

fineness = 0.1

variant_ <- unique(Korea_dominance$classified_label)[1]

N <- df_variant_$numerator

date_ <- df_variant_$date

df_variant_ <- Korea_data %>% 
  dplyr::filter(classified_label == variant_)

# break down the function
date <- as_date(date_)
Date_numeric <- as.numeric(date - min(date))

d <- as.factor(weekdays(date))
this_k <-  floor(length(Date_numeric) / denominator)

the_formula <- switch(dow,
                      none = as.formula(N ~ s(Date_numeric, k = this_k, bs = bs)),
                      re = as.formula(N ~ s(Date_numeric, k = this_k, bs = bs) + s(d, bs = 're')),
                      fixed = as.formula(N ~ s(Date_numeric, k = this_k, bs = bs) + d))

this_model <- gam(formula = the_formula,
                  family = 'nb')

data.frame(Date_numeric = seq(min(Date_numeric), max(Date_numeric), fineness)) %>%
  dplyr::mutate(date = Date_numeric + min(date),
         d = as.factor(weekdays(date))) -> new_data0

full_predictions <- predict(object = this_model,
                            type = 'link',
                            se.fit = TRUE,
                            newdata = new_data0,
                            exclude = 's(d)')

new_data0 %>%
  dplyr::mutate(full_fit = exp(full_predictions$fit),
         full_upper = exp(full_predictions$fit + 2 * full_predictions$se.fit),
         full_lower = exp(full_predictions$fit - 2 * full_predictions$se.fit))  -> model_fit

# compare the fitness

raw_ <- data.frame(date = date, N = N, type = "raw")
fit_ <- data.frame(date = model_fit$date, N = model_fit$full_fit, type = "fit")


ggplot() + 
  geom_point(data = raw_, mapping = aes(x = date, y = N), color = "blue4", alpha = .5) + 
  geom_line(data = fit_, mapping = aes(x = date, y = N))


################################################################################
matrix0 <- predict(object = this_model,
                   newdata = new_data0,
                   type = 'lpmatrix',
                   exclude = 's(d)')


new_data0 %>%
  dplyr::mutate(Date_numeric = Date_numeric + 1e-7) -> new_data1

matrix1 <- predict(object = this_model,
                   newdata = new_data1,
                   type = 'lpmatrix',
                   exclude = 's(d)')

X <- (matrix1 - matrix0) / 1e-7

new_data0 %>%
  dplyr::mutate(rate = as.vector(X %*% coef(this_model)),
         upper_rate = as.vector(rate + 2 * rowSums(X %*% this_model$Vp * X) ^ 0.5),
         lower_rate = as.vector(rate - 2 * rowSums(X %*% this_model$Vp * X) ^ 0.5)) -> fitted_rates


ggplot(fitted_rates, aes(x = date, y = rate)) + geom_line()


