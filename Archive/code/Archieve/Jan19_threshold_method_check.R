# Load required libraries
library(dplyr)
library(glasso)
library(genlasso)
library(lubridate)
library(signal)
library(TTR)

# Source required custom functions
source('Background/functions_0706.R')
source('Code/functions_franky_0725.R')
source('Background/growth_rate_function_original.R')

# Load dataset
url <- "UKHSA-UConn-variant-modelling/variant_modelling/data/summary_GISAID_20240918.csv"
data <- read.csv(url)

# Extract top 30 countries by total infections
main_countries <- data %>%
  group_by(country) %>%
  summarise(num_total = sum(numerator)) %>%
  arrange(desc(num_total)) %>%
  head(30) %>%
  pull(country)

# Filter data for main countries
data2 <- data %>%
  dplyr::filter(country %in% main_countries)

# Count unique lineages in the filtered dataset
num_lineages <- length(unique(data2$lineage))
print(paste("Number of unique lineages:", num_lineages)) # Output: 3,746

# Load alias list for lineage classification
alias_list <- get_alias()

# Rename columns and select relevant variables
data3 <- data2 %>%
  rename(n = numerator) %>%
  select(country, date, lineage, n)


save_url <- "Code/ProcessedData/global_classification_to_append.rds"

classification_to_append <- readRDS(save_url)

classified_data <- data3 |>
  dplyr::mutate(
    decimal_lineage = paste0(lineage, '.')
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
    length_lineage = stringr::str_count(lineage, '\\.'),
    classified_label = dplyr::if_else((length_lineage - length_class) <= 1,
                                      classified_label,
                                      NA),
    classified_label = dplyr::if_else(is.na(classified_label),
                                      "Other",
                                      classified_label)
  ) |>
  dplyr::filter(classified_label != 'Other') |>
  dplyr::select(country, date, lineage, n, classified_label)


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

classified_data$date <- as.Date(classified_data$date)

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

classified_data$share[is.nan(classified_data$share)] <- 0

# do the histogram plots 
library(dplyr)
library(ggplot2)
library(RColorBrewer)

create_country_plot <- function(classified_data, select) {
  
  # 生成颜色对应表，确保颜色一致
  color_palette <- c(
    brewer.pal(9, "Set1"), 
    brewer.pal(8, "Dark2"), 
    brewer.pal(5, "Set3")[1:5]
  )
  unique_labels <- unique(classified_data$classified_label)
  label_colors <- setNames(color_palette[1:length(unique_labels)], unique_labels)
  # 数据处理
  processed_data <- classified_data %>%
    dplyr::filter(country == select) %>%
    dplyr::mutate(year_month = format(as.Date(date), "%Y-%m")) %>%
    dplyr::group_by(year_month, classified_label) %>%
    dplyr::summarize(
      total_numerator = sum(numerator, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(year_month) %>%
    dplyr::mutate(
      total_all_numerator = sum(total_numerator, na.rm = TRUE),
      share = total_numerator / total_all_numerator
    ) %>%
    dplyr::ungroup()
  
  # 生成堆叠柱状图
  ggplot(processed_data, aes(x = year_month, y = share, fill = classified_label)) +
    geom_bar(stat = "identity", position = "stack", color = "black", size = 0.25) +
    scale_fill_manual(values = label_colors) + # 确保颜色一致
    labs(
      title = paste("Monthly Share of", select),
      x = "",
      y = "Share"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom",
      legend.title = element_blank(), 
      legend.text = element_text(size = 8),
      legend.key.width = unit(1.5, "lines"), 
      legend.spacing.x = unit(0.3, "cm"),  
      legend.box = "horizontal",  
      legend.box.just = "center",
      panel.grid = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 3)) + 
    scale_x_discrete(
      breaks = processed_data$year_month[grep("-01$|-07$", processed_data$year_month)], 
      labels = function(x) format(as.Date(paste0(x, "-01")), "%Y-%b")
    )
}



create_country_plot(classified_data, select = "South Korea")
create_country_plot(classified_data, select = "United Kingdom")
create_country_plot(classified_data, select = "USA")
create_country_plot(classified_data, select = "Australia")
create_country_plot(classified_data, select = "Brazil")
create_country_plot(classified_data, select = "Canada")























classified_label.select <- classified_label.appearance %>% 
  dplyr::filter(count > 130 & count_numerator_gt_1 >= 60) %>% 
  ungroup() %>% 
  mutate(country_label = paste0(country, sep = "-", classified_label)) %>% 
  dplyr::select(country_label)

classified_data <- classified_data %>% 
  mutate(country_label =  paste0(country, sep = "-", classified_label)) %>% 
  dplyr::filter(country_label %in% classified_label.select$country_label) %>% 
  dplyr::select(country_label, date, numerator, denominator, share)


# initial stage after adjustment
filtered_data <- classified_data %>%
  arrange(country_label, date) %>% 
  group_by(country_label) %>%
  mutate(
    numerator_positive = numerator > 0,
    run_length = with(rle(numerator_positive), rep(lengths, lengths)),
    valid_period = numerator_positive & run_length >= 7
  ) %>%
  dplyr::filter(cumsum(valid_period) > 0) %>%
  dplyr::select(-numerator_positive, -run_length, -valid_period) %>%
  ungroup()

data_list <- list()

for (i in unique(filtered_data$country_label)) {
  
  # i = "Russia-B.1.1" # for test
  variant_data <- filtered_data |> 
    dplyr::filter(country_label == i) |> 
    arrange(country_label, date) |> 
    mutate(
      sharing = dplyr::coalesce(share, 0),
      gap = as.numeric(date - dplyr::first(date))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-share)
  
  data_list[[i]] <- variant_data
  
}


################################################################################
# start to define the thresholds
################################################################################

count_gaps_over_threshold <- function(data_list, share_threshold) {
  
  # share_threshold = .4
  
  results_df <- data.frame(var = character(), 
                           date = as.Date(character()), 
                           gap = numeric(), 
                           total_gap = numeric(), 
                           stringsAsFactors = FALSE)
  
  for (var in names(data_list)) {
    df <- data_list[[var]]
    
    result <- df %>%
      dplyr::filter(sharing >= share_threshold) %>%
      slice(1) %>%
      dplyr::select(date, gap)
    
    if (nrow(result) == 0) {
      result <- data.frame(date = NA, gap = NA)
    }
    
    result$var <- var
    result$total_gap <- tail(df$gap, 1)
    results_df <- rbind(results_df, result)
  }
  
  return(results_df)
}

result <- count_gaps_over_threshold(data_list, share_threshold = 0.5)

sum(!is.na(result$gap))
sum(result$gap > 14, na.rm = TRUE)











