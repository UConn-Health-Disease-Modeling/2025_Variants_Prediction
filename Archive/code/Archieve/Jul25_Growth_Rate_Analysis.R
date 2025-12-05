# set root directory
setwd("/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis")

# source the functions (functions_0706.R)
source('Background/functions_0706.R')

# source my functions
source('Code/functions_franky_0725.R')

# load packages 
library(magrittr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(tidyr)
library(purrr)
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(zoo)
library(wavelets)
library(nlme)
library(car)
# install.packages("psych")
library(psych)

# load the data
output_url <- "Data/output_SouthKorea2.rds"
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

# calculate the accumulated monthly share
dominant_index <- monthly_output %>% 
  dplyr::group_by(classified_label) %>% 
  dplyr::summarise(
    acc_sum = sum(share, na.rm = TRUE), 
    .groups = 'drop'
  ) %>% 
  dplyr::arrange(
    desc(acc_sum)
  ) 

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
  )

####################################
######## 1st: Based on Time ########
####################################

# extract first 30 days data
output_90d <- output %>%
  dplyr::filter(
    days_since_first <= 30
  ) %>%
  dplyr::distinct() %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::slice(-1) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    classified_label, gr, days_since_first
  )

output_90d_wide <- output_90d %>%
  pivot_wider(
    names_from = days_since_first,
    values_from = gr
  ) %>% 
  left_join(
    dominant_index, 
    by = "classified_label"
  )

# convert the data to a matrix for convenience 
data <- as.matrix(output_90d_wide[, 2:32])
target <- as.matrix(output_90d_wide[, 33])

# compress the data (91 -> 9)
compressed_data <- t(apply(data, 1, paa_transform, num_segments = 9))

# take a example to compare the original and compressed data
index <- sample(1:29, 1)
selected_row <- data[index, ]
compressed_row <- compressed_data[index, ]

# prepare the data
original_df <- data.frame(
  Time = 1:length(selected_row),
  Value = selected_row,
  Type = "Original"
)

compressed_df <- data.frame(
  Time = seq(1, length(selected_row), length.out = length(compressed_row)),
  Value = compressed_row,
  Type = "Compressed"
)

# draw the plot
ggplot() +
  geom_line(data = original_df, aes(x = Time, y = Value, color = Type), size = 1) +
  geom_line(data = compressed_df, aes(x = Time, y = Value, color = Type), size = 1, linetype = "dashed") +
  labs(title = "Comparison of Original and Compressed Time Series",
       x = "Time",
       y = "Value") +
  scale_color_manual(values = c("Original" = "blue", "Compressed" = "red")) +
  theme_minimal()

## Do the Linear regression model 
target <- target[1:29]

compressed_data_df <- as.data.frame(compressed_data)
compressed_data_df$target <- target

model <- lm(target ~ ., data = compressed_data_df)
summary(model)

results <- correlation_tests(compressed_data_df, "target")
print(results)


####################################
##### 2nd: Based on Prevalence #####
####################################
output_filtered <- output %>%
  dplyr::group_by(
    classified_label
  ) %>%
  dplyr::filter(
    any(prevalence >= 0.05)
  )


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

table(output_initial$classified_label)

# hypo1: the time to the threshold (0.05)
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


# hypo2: average growth rate at initial stage 
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


output_initial <- output_initial %>% 
  select(
    classified_label, days_since_first, gr
  )

output_paa <- output_initial %>%
  group_by(classified_label) %>%
  summarise(paa_gr = list(PAA(gr, 7))) %>%
  unnest(cols = c(paa_gr)) %>% 
  mutate(days = rep(1:7, 23))


output_paa_wide <- output_paa %>%
  pivot_wider(
    names_from = days,
    values_from = paa_gr
  ) %>% 
  left_join(
    dominant_index, 
    by = "classified_label"
  )

model2 <- lm(acc_sum ~ ., 
             data = output_paa_wide %>% 
               select(-classified_label))
summary(model2)

# consider the auto-correlation problem 
# install.packages("keras")
library(keras)

# Install TensorFlow if not already installed
install_keras()

timesteps <- 50 

# Define LSTM model
model <- keras_model_sequential() %>%
  layer_lstm(units = 50, input_shape = c(timesteps, 1), return_sequences = FALSE) %>%
  layer_dense(units = 1)

model %>% compile(
  loss = 'mean_squared_error',
  optimizer = 'adam'
)

input_matrix <- as.matrix(output_paa_wide[, 2:8])
output_matrix <- as.matrix(output_paa_wide[, 9])
  
history <- model %>% fit(
  input_matrix, output_matrix,
  epochs = 50,
  batch_size = 8,
  validation_split = 0.2
)



####################################
####### 3rd: Extend the Time #######
####################################

# first 30 days data
model_30days <- perform_pca_regression(data = output, days_to_use = 30, freq = 5)
model_45days <- perform_pca_regression(data = output, days_to_use = 45, freq = 5)
model_60days <- perform_pca_regression(data = output, days_to_use = 60, freq = 5)


summary(model_30days)
summary(model_45days)
summary(model_60days)

