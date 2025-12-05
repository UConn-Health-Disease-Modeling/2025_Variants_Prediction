# source the functions (functions_0706.R)
source(glue::glue('/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis/Background/functions_0706.R'))
source('Code/functions_franky_0725.R')

################################################################################
# DON'T RUN
# list.files()
url <- "gisaid_variants_statistics_2024_08_05_1942/gisaid_variants_statistics.json"
data <- jsonlite::fromJSON(url)$stats

data_df <- data.frame(
  count = numeric(),
  value = character(),
  date = as.Date(character()),
  country = character()
)

for (date in names(data)) {
  
  # names(data) # do not run
  # date <- "2020-10-18"  # for test
  
  cat("processing", date, "\n")
  
  list_of_date <- data[[date]]
  
  for (country in names(list_of_date)) {
    
    # country <- "United Kingdom" # for test
    # cat(country, "\n")
    
    country_of_date <- list_of_date[[country]]
    
    to_append <- country_of_date$submissions_per_lineage
    
    if(length(to_append) > 0){
      
      to_append$date <- date
      to_append$country <- country
      
      data_df <- rbind(data_df, to_append)
      
    } 
    
  }
  
}

colnames(data_df) <- c("count", "lineage", "date", "country")

# save the dataset
save_url <- "Data/gisaid_variants_statistics.rds"
saveRDS(data_df, save_url)
################################################################################



library(lubridate)
df <- readRDS("Data/gisaid_variants_statistics.rds")

length(unique((df |> filter(country == "United Kingdom"))$lineage))



df$date <- as.Date(df$date, format = "%Y-%m-%d")

# filter existing lineages (6 months period from latest date Aug 04, 2024)
lineage_select <- df |>
  dplyr::group_by(lineage) |>
  dplyr::summarise(
    last_appreance = max(date),
    .groups = "drop"
  ) |>
  dplyr::filter(last_appreance <= "2024-02-04") |>
  dplyr::pull(lineage)

df1 <- df |>
  dplyr::filter(
    lineage %in% lineage_select
  ) |>
  dplyr::group_by(
    country, 
    lineage
  ) |> 
  tidyr::complete(
    date = seq(min(date), max(date), by = "7 days")
  ) |> 
  tidyr::replace_na(
    list(count = 0)
  ) |> dplyr::ungroup() |> 
  dplyr::group_by(
    country, date
  ) |>
  dplyr::mutate(
    total_count_at_date = sum(count, na.rm = TRUE)
  ) |>
  dplyr::mutate(
    share = count / total_count_at_date
  ) |> dplyr::select(country, date, lineage, count, total_count_at_date, share) |> 
  dplyr::mutate(
    share = ifelse(is.nan(share), 0, share)
  )


# "China" %in% unique(df1$country)

# select the main countries 
df2 <- df1 |> 
  dplyr::group_by(
    country
  ) |> 
  dplyr::mutate(
    sum_count = sum(count, na.rm = TRUE)
  ) |>
  dplyr::ungroup() |> 
  dplyr::filter(sum_count >= 3000) |>
  # dplyr::filter(country %in% c("United Kingdom", "South Korea", "USA", "India", "Japan", "Canada", 
  #                              "Australia", "France", "Mexico", "Italy", "Spain")) |> 
  dplyr::arrange(country, lineage, date) |>
  dplyr::select(
    country, date, lineage, count, total_count_at_date, share
  )

length(unique(df2$country))



length(unique(df2$country)) # keep 64 countries
length(unique(df2$lineage)) # 605 lineages

# select the lineage with duration >= 60 days 
df2 |> dplyr::group_by(
  country, lineage
  ) |> 
  dplyr::summarise(
    start_date = min(date, na.rm = TRUE),
    end_date = max(date, na.rm = TRUE),
    duration = as.numeric(difftime(end_date, start_date, units = "days")),
    .groups = "drop"
  ) |> 
  dplyr::filter(
    duration >= 60
  ) -> lineage_select

df3 <- df2 |>
  dplyr::left_join(
    lineage_select, 
    by = c("country", "lineage")
  ) |> 
  dplyr::filter(
    !is.na(duration)
  )

# select the lineage with at least 3 non-zero values during first 6 weeks
df3 |>
  dplyr::group_by(
    country, lineage
  ) |>
  dplyr::arrange(
    date
  ) |>  
  dplyr::slice_head(
    n = 8
  ) |> 
  dplyr::summarise(
    non_zero_count = sum(count > 0, na.rm = TRUE), .groups = "drop"
  ) |>  
  dplyr::filter(
    non_zero_count > 4
  ) -> lineage_select

df4 <- df3 |>
  dplyr::left_join(
    lineage_select, 
    by = c("country", "lineage")
  ) |> 
  dplyr::filter(
    !is.na(non_zero_count)
  )

# # select the lineage with at least 3 non-zero values during first 6 weeks
# df4 |> 
#   dplyr::group_by(
#     country, lineage
#   ) |> 
#   dplyr::summarise(
#     max_share = max(share), 
#     .groups = "drop"
#   ) |> 
#   dplyr::filter(
#     max_share >= .1
#   ) -> lineage_select
# 
# 
# df5 <- df4 |>
#   dplyr::left_join(
#     lineage_select, 
#     by = c("country", "lineage")
#   ) |> 
#   dplyr::filter(
#     !is.na(max_share)
#   )

# calculate the accumulative share
df4 |> 
  dplyr::group_by(
    country, lineage
  ) |>
  dplyr::summarise(
    acc_share = sum(share), 
    .groups = "drop"
  ) -> response_df

df4 |> 
  dplyr::group_by(
    country, lineage
  ) |>
  dplyr::arrange(
    date
  ) |>  
  dplyr::slice_head(
    n = 8
  ) -> input_df_8

df4 |> 
  dplyr::group_by(
    country, lineage
  ) |>
  dplyr::arrange(
    date
  ) |>  
  dplyr::slice_head(
    n = 6
  ) -> input_df_6

df4 |> 
  dplyr::group_by(
    country, lineage
  ) |>
  dplyr::arrange(
    date
  ) |>  
  dplyr::slice_head(
    n = 4
  ) -> input_df_4
  
library(dplyr)
library(tidyr)

input_df_8 <- input_df_8 %>%
  group_by(country, lineage) %>%
  mutate(id = row_number()) %>%
  ungroup()

count_8 <- input_df_8 %>%
  select(country, lineage, id, count) %>%
  pivot_wider(
    names_from = id,
    values_from = count,
    names_prefix = "count_"
  )

dataset_8weeks <- input_df_8 %>%
  select(country, lineage, id, share) %>%
  pivot_wider(
    names_from = id,
    values_from = share,
    names_prefix = "share_"
  ) %>% 
  left_join(count_8, by = c("country", "lineage")) %>%
  left_join(response_df, by = c("country", "lineage"))



input_df_6 <- input_df_6 %>%
  group_by(country, lineage) %>%
  mutate(id = row_number()) %>%
  ungroup()

count_6 <- input_df_6 %>%
  select(country, lineage, id, count) %>%
  pivot_wider(
    names_from = id,
    values_from = count,
    names_prefix = "count_"
  )

dataset_6weeks <- input_df_6 %>%
  select(country, lineage, id, share) %>%
  pivot_wider(
    names_from = id,
    values_from = share,
    names_prefix = "share_"
  ) %>% 
  left_join(count_6, by = c("country", "lineage")) %>%
  left_join(response_df, by = c("country", "lineage"))



input_df_4 <- input_df_4 %>%
  group_by(country, lineage) %>%
  mutate(id = row_number()) %>%
  ungroup()

count_4 <- input_df_4 %>%
  select(country, lineage, id, count) %>%
  pivot_wider(
    names_from = id,
    values_from = count,
    names_prefix = "count_"
  )

dataset_4weeks <- input_df_4 %>%
  select(country, lineage, id, share) %>%
  pivot_wider(
    names_from = id,
    values_from = share,
    names_prefix = "share_"
  ) %>% 
  left_join(count_4, by = c("country", "lineage")) %>%
  left_join(response_df, by = c("country", "lineage"))

# save the dataset
save_url <- "Data/machine_learning_datasets.rds"
saveRDS(list(four = dataset_4weeks, 
             six = dataset_6weeks, 
             eight = dataset_8weeks), 
        save_url)

datasets <- readRDS(save_url)


t1 <- datasets$four |> group_by(country) |>
  mutate(mean_share = (share_1 + share_2 + share_3 + share_4)/4, 
         cor4 = cor(acc_share, mean_share)) |> 
  select(country, cor4) |> unique()

t2 <- datasets$six |> group_by(country) |>
  mutate(mean_share = (share_1 + share_2 + share_3 + share_4 + share_5 + share_6)/6, 
         cor6 = cor(acc_share, mean_share)) |> 
  select(country, cor6) |> unique()


t3 <- datasets$eight |> group_by(country) |>
  mutate(mean_share = (share_1 + share_2 + share_3 + share_4 + share_5 + share_6 + share_7 + share_8)/8, 
         cor8 = cor(acc_share, mean_share)) |> 
  select(country, cor8) |> unique()

summary <- t1 |> 
  left_join(t2, by = "country") |> 
  left_join(t3, by = "country")

top_15_countries <- c(
  "United States", "China", "India", "Japan", "Germany", 
  "United Kingdom", "France", "Brazil", "Italy", "Canada", 
  "Russia", "Australia", "South Korea", "Saudi Arabia", "Mexico"
)


summary_top15 <- summary |> 
  filter(country %in% top_15_countries)

datasets1 <- datasets$four |> filter(country %in% top_15_countries) |> 
  select(country, share_1, share_2, share_3, share_4, acc_share)
datasets2 <- datasets$six |> filter(country %in% top_15_countries) |> 
  select(country, share_1, share_2, share_3, share_4, share_5, share_6, acc_share)
datasets3 <- datasets$eight |> filter(country %in% top_15_countries) |> 
  select(country, share_1, share_2, share_3, share_4, share_5, share_6, share_7, share_8, acc_share)



r_squared_results <- data.frame(country = unique(datasets1$country), 
                                r_squared_4 = NA, 
                                r_squared_6 = NA, 
                                r_squared_8 = NA)

for (i in seq_along(r_squared_results$country)) {
  
  country_data1 <- datasets1 %>% filter(country == r_squared_results$country[i])
  country_data2 <- datasets2 %>% filter(country == r_squared_results$country[i])
  country_data3 <- datasets3 %>% filter(country == r_squared_results$country[i])
  
  model1 <- lm(acc_share ~ share_1 + share_2 + share_3 + share_4, 
               data = country_data1)
  model2 <- lm(acc_share ~ share_1 + share_2 + share_3 + share_4 + share_5 + share_6, 
               data = country_data2)
  model3 <- lm(acc_share ~ share_1 + share_2 + share_3 + share_4 + share_5 + share_6 + share_7 + share_8, 
               data = country_data3)
  
  r_squared_results$r_squared_4[i] <- summary(model1)$r.squared
  r_squared_results$r_squared_6[i] <- summary(model2)$r.squared
  r_squared_results$r_squared_8[i] <- summary(model3)$r.squared
}

write.csv(r_squared_results, file = "r_squared_results_lm.csv", row.names = FALSE)

library(e1071) 

r_squared_results_svm <- data.frame(country = unique(datasets1$country), 
                                    r_squared_4 = NA, 
                                    r_squared_6 = NA, 
                                    r_squared_8 = NA)

# Define a function to calculate R-squared for SVM
calculate_r_squared <- function(true_values, predicted_values) {
  rss <- sum((true_values - predicted_values)^2)
  tss <- sum((true_values - mean(true_values))^2)
  return(1 - rss/tss)
}

# Loop through each country and fit the SVM model
for (i in seq_along(r_squared_results_svm$country)) {
  
  country_data1 <- datasets1 %>% filter(country == r_squared_results_svm$country[i])
  country_data2 <- datasets2 %>% filter(country == r_squared_results_svm$country[i])
  country_data3 <- datasets3 %>% filter(country == r_squared_results_svm$country[i])
  
  # Fit the SVM model with 4 predictors
  svm_model1 <- svm(acc_share ~ share_1 + share_2 + share_3 + share_4, 
                    data = country_data1)
  predicted1 <- predict(svm_model1, country_data1)
  
  # Fit the SVM model with 6 predictors
  svm_model2 <- svm(acc_share ~ share_1 + share_2 + share_3 + share_4 + share_5 + share_6, 
                    data = country_data2)
  predicted2 <- predict(svm_model2, country_data2)
  
  # Fit the SVM model with 8 predictors
  svm_model3 <- svm(acc_share ~ share_1 + share_2 + share_3 + share_4 + share_5 + share_6 + share_7 + share_8, 
                    data = country_data3)
  predicted3 <- predict(svm_model3, country_data3)
  
  # Calculate R-squared for each model
  r_squared_results_svm$r_squared_4[i] <- calculate_r_squared(country_data1$acc_share, predicted1)
  r_squared_results_svm$r_squared_6[i] <- calculate_r_squared(country_data2$acc_share, predicted2)
  r_squared_results_svm$r_squared_8[i] <- calculate_r_squared(country_data3$acc_share, predicted3)
}


write.csv(r_squared_results_svm, file = "r_squared_results_svm.csv", row.names = FALSE)
