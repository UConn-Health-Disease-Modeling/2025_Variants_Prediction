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
# p_lim <- 0.01
# 
# # set n limit
# n_lim <- 100
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
# saveRDS(classification_to_append, "Data/Aug15_classification_to_append_Combined.rds")

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

rm(alias_list, classified_data, data, append_url, url)





################################################################################
# 30 Days data
UK_data_30_ <- UK_data %>% 
  dplyr::filter(
    gap <= 30
  ) %>% 
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    sum(sharing == 0) <= (30-3)
  ) %>%
  dplyr::ungroup()

name_to_keep_UK_ <- unique(UK_data_30_$classified_label)

Korea_data_30_ <- Korea_data %>% 
  dplyr::filter(
    gap <= 30
  ) %>% 
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    sum(sharing == 0) <= (30-3)
  ) %>%
  dplyr::ungroup()

name_to_keep_Korea_ <- unique(Korea_data_30_$classified_label)





################################################################################
# apply to all the data 
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
    
    # Process 'numerator' and 'sharing' columns
    abs_data <- process_growth_rates(variant_df_, "numerator")
    share_data <- process_growth_rates(variant_df_, "sharing")
    
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





Korea_list_30_ <- generate_growth_rate_list(Korea_data_30_)
UK_list_30_ <- generate_growth_rate_list(UK_data_30_)



# 60 Days data
UK_data_60_ <- UK_data %>% 
  dplyr::filter(
    gap <= 60
  ) %>% 
  dplyr::filter(
    classified_label %in% name_to_keep_UK_
  )


Korea_data_60_ <- Korea_data %>% 
  dplyr::filter(
    gap <= 60
  ) %>% 
  dplyr::filter(
    classified_label %in% name_to_keep_Korea_
  )

Korea_list_60_ <- generate_growth_rate_list(Korea_data_60_)
UK_list_60_ <- generate_growth_rate_list(UK_data_60_)

r_name_ <- sample(unique(Korea_data_30_$classified_label), 1)
# r_name_ <- "HK.3 (XBB.1.9.2.5.1.1.3)"

# r_raw <- Korea_list_30_[[r_name_]]$raw
# r_fit <- Korea_list_30_[[r_name_]]$fit
r_raw <- Korea_list_60_[[r_name_]]$raw
r_fit <- Korea_list_60_[[r_name_]]$fit

p1 <- ggplot() +
  geom_line(data = r_fit, mapping = aes(x = date, y = sharing_fit)) +
  geom_point(data = r_raw, mapping = aes(x = date, y = sharing), color = "red4") + 
  theme_bw() +
  ggtitle(bquote(bold("Share Fit: ") ~ .(r_name_))) + xlab("Date") + ylab("Share")
  theme(plot.title = element_text(face = "bold"))

p2 <- ggplot() +
  geom_line(data = r_fit, mapping = aes(x = date, y = numerator_fit)) +
  geom_point(data = r_raw, mapping = aes(x = date, y = numerator), color = "blue4") + 
  theme_bw() +
  ggtitle(bquote(bold("Absolute Fit: ") ~ .(r_name_))) + xlab("Date") + ylab("Number of Infection")
  theme(plot.title = element_text(face = "bold"))

gridExtra::grid.arrange(p1, p2, nrow = 1)

Korea_dominance %>% 
  dplyr::filter(
    classified_label %in% name_to_keep_Korea_
  ) %>% 
  dplyr::mutate(
    dominance_indice = ifelse(dominance > 40, "yes", "no")
  ) -> Korea_dominance

dom_Korea_names <- Korea_dominance$classified_label[Korea_dominance$dominance_indice == "yes"]
non_dom_Korea_names <- Korea_dominance$classified_label[Korea_dominance$dominance_indice == "no"]


list_30_ <- list(
  Korea = Korea_list_30_ , 
  UK = UK_list_30_
)

list_60_ <- list(
  Korea = Korea_list_60_ , 
  UK = UK_list_60_
)

list_dom <- list(
  Korea = Korea_dominance %>% dplyr::filter(classified_label %in% name_to_keep_Korea_), 
  UK = UK_dominance %>% dplyr::filter(classified_label %in% name_to_keep_UK_)
)

# # check the current root 
# list.files()
# saveRDS(list_30_, "Code/ProcessedData/list_30_.rds")
# saveRDS(list_60_, "Code/ProcessedData/list_60_.rds")
# saveRDS(list_dom, "Code/ProcessedData/list_dom_.rds")



