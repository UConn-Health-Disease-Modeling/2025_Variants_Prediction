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

# Summarise total infections by country and lineage
summarised <- data3 %>%
  group_by(country, lineage) %>%
  summarise(total_infections = sum(n), .groups = "drop")

# #
# # set p limit
# p_lim <- 0.01
# #
# # set n limit
# n_lim <- 100
# #
# # # get the classification to append
# # # save the result since the function costs a lot of time
# classification_to_append <- get_classification(
#   lineages = summarised$lineage,
#   number_sequences = summarised$n,
#   p_lim = p_lim, n_lim = n_lim, alias_list
# ) |>
#   dplyr::mutate(
#     decimal_lineage = paste0(classified_unasliased, '.')
#   )



save_url <- "Code/ProcessedData/global_classification_to_append.rds"
# saveRDS(classification_to_append, save_url)


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

# table(data3$country)
# Australia        Austria        Belgium         Brazil         Canada          China Czech Republic        Denmark 
# 38875          21612          27369          21508          71208          17429          12102          41173 
# France        Germany          India        Ireland         Israel          Italy          Japan     Luxembourg 
# 46806          68011          20303          20012          21134          32102          57363          15261 
# Mexico    Netherlands         Norway           Peru         Poland         Russia       Slovenia    South Korea 
# 17947          34404          15340          12620          13100          15170          14583          33927 
# Spain         Sweden    Switzerland         Turkey United Kingdom            USA 
# 42678          36144          26652           5776          99752         185312 

# --->

# table(classified_data$country)
# Australia        Austria        Belgium         Brazil         Canada          China Czech Republic        Denmark 
#     12499           8448          12865           9762          24821           4517           6831          17079 
#    France        Germany          India        Ireland         Israel          Italy          Japan     Luxembourg 
#     17930          28132          10075           8534          10179          13593          16707           7491 
#    Mexico    Netherlands         Norway           Peru         Poland         Russia       Slovenia    South Korea 
#      9618          14879           8657           5349           6935           8134           7426           9941 
#     Spain         Sweden    Switzerland         Turkey United Kingdom            USA 
#     16245          14579          14514           3222          42522          72354 

# classified_data |> dplyr::filter(country == "Australia")

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

classified_label.appearance <- classified_data %>% 
  group_by(country, classified_label) %>% 
  summarise(
    count = n(),
    count_numerator_gt_1 = sum(numerator > 1)
  )




# the 5% position is 130, then dplyr::filter all those appearance less than 130
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

# pdf("Manuscripts/plots/initial_stage_adjustment.pdf", width = 10, height = 8) # Specify file name and size
# 
# # Set plotting layout for 4x2 grid
# par(mfrow = c(4, 2))
# 
# # Set seed for reproducibility
# set.seed(27)
# 
# # Sample two labels from South Korea
# label.Korea <- classified_data %>%
#   filter(startsWith(country_label, "South Korea")) %>%
#   sample_n(2) %>%
#   pull(country_label)
# 
# # Sample two labels from the United Kingdom
# label.UK <- classified_data %>%
#   filter(startsWith(country_label, "United Kingdom")) %>%
#   sample_n(2) %>%
#   pull(country_label)
# 
# # Function to plot share data for a specific country label
# plot_share <- function(data, label, xlab = "Time (Day)", ylab = "Share (Country-level)") {
#   # Remove parentheses unless the content is in the format B.x.x.x
#   label_clean <- ifelse(grepl("\\(B\\.[0-9]+(\\.[0-9]+)*\\)", label),
#                         label,
#                         gsub("\\s*\\(.*\\)", "", label))
#   
#   # Extract the data for the specified label
#   share_data <- data %>% filter(country_label == label) %>% pull(share)
#   
#   # Create the plot
#   plot(
#     share_data,
#     type = "p", pch = 3, cex = 0.6, xlab = xlab, ylab = ylab,
#     main = bquote(bold(.(label_clean))) # Add bold title with cleaned label
#   )
#   
#   # Add gray rectangle for days 0-60
#   rect(xleft = 0, xright = 30, ybottom = min(share_data, na.rm = TRUE), 
#        ytop = max(share_data, na.rm = TRUE), col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
# }
# 
# 
# 
# # Plot data for selected labels
# for (label in c(label.Korea, label.UK)) {
#   plot_share(classified_data, label)
#   plot_share(filtered_data, label)
# }
# 
# # Close the PDF device
# dev.off()



# still have 511 variants among the 30 countries 
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

data_list.days <- list()

for (j in 30:60) {
  
  Index <- paste0("days", sep = "_", j)
  
  data_list.days[[Index]] <- list()
  
  for (i in names(data_list)) {
    
    data_list.days[[Index]][[i]] <- data_list[[i]] %>% 
      dplyr::filter(gap <= j) 
    
  }

}

lable.eg <- sample(names(data_list.days$days_6), 1)
plot(data_list.days[["days_60"]][[lable.eg]] %>% dplyr::select(numerator) %>% unlist(), 
     type = "p", pch = 3, cex = .6, xlab = "days", ylab = "frequences")

rm(list = setdiff(ls(), ls(pattern = "^data")))



# data_list.days.smoothed <- data_list.days
# 
# for (i in names(data_list.days)) {  
#   
#   # i <- "days_60"
#   cat("processing", i, "\n")
#   
#   data <- data_list.days[[i]]
#   
#   for (j in names(data)) {
#     
#     # j = "Poland-BA.2 (BA.2)"
#     
#     # cat("processing", i, "---", j, "\n")
#     
#     knots <- as.numeric(sub("days_", "", i)) %>% sqrt() %>% round()
#     
#     sub_data <- data[[j]]
#     
#     # Extrapolation the sequence to reduce the marginal effects 
#     window_size <- 5
#     moving_avg_sharing <- SMA(sub_data$sharing, n = window_size)
#     moving_avg_freq    <- SMA(sub_data$numerator, n = window_size)
#     
#     last_sharing_avg <- tail(moving_avg_sharing, 5)
#     add_sharing <- rep(last_sharing_avg, round(dim(sub_data)[1]/10))
#     last_freq_avg <- tail(moving_avg_freq, 5)
#     add_freq <- rep(last_freq_avg, round(dim(sub_data)[1]/10))
#     
#     extend_sharing <- c(sub_data$sharing, add_sharing)
#     extend_freq    <- c(sub_data$numerator, add_freq)
#       
#     smoothed_sharing <- sgolayfilt(extend_sharing, p = 3, n = 45)
#     smoothed_freq    <- sgolayfilt(extend_freq, p = 3, n = 45)
#     
#     new.grid <- data.frame(gap = seq(min(sub_data$gap), max(sub_data$gap + length(add_sharing)), length.out = 100*(length(extend_sharing)[1] - 1)))
#     new.grid$country_label <- j
#     
#     org_gap <- seq(from = 0, to = length(extend_sharing) - 1)
#     
#     new.grid$smoothed_sharing <- approx(x = org_gap, y = smoothed_sharing, xout = new.grid$gap)$y
#     new.grid$smoothed_freq    <- approx(x = org_gap, y = smoothed_freq, xout = new.grid$gap)$y
#     
#     data_list.days.smoothed[[i]][[j]] <- new.grid[1:(100*(dim(sub_data)[1]-1)), ]
#     
#   }
#   
# }


# list.files()
S_G.smoothed.list.url <- "Code/ProcessedData/S_G.smoothed.list.rds"
# saveRDS(data_list.days.smoothed, S_G.smoothed.list.url)


data_list.days.smoothed <- readRDS(S_G.smoothed.list.url)

# smoothed.gr <- data_list.days.smoothed
# loess_span <- 0.3  
# 
# for (i in names(data_list.days.smoothed)) {
#   
#   cat("processing-", i, "\n")
#   
#   i = "days_60"
#   
#   for (j in names(data_list.days.smoothed[[i]])) {
#     
#     df <- data_list.days.smoothed[[i]][[j]]
#     
#     smoothed_sharing_growth_raw = c(NA, diff(df$smoothed_sharing) / lag(df$smoothed_sharing)[-1])
#     smoothed_freq_growth_raw = c(NA, diff(df$smoothed_freq) / lag(df$smoothed_freq)[-1])
#     
#     loess_sharing_model <- loess(smoothed_sharing_growth_raw ~ seq_along(smoothed_sharing_growth_raw), span = loess_span, na.action = na.exclude)
#     loess_freq_model <- loess(smoothed_freq_growth_raw ~ seq_along(smoothed_freq_growth_raw), span = loess_span, na.action = na.exclude)
#     
#     smoothed.gr[[i]][[j]]$smoothed_sharing_growth <- predict(loess_sharing_model)
#     smoothed.gr[[i]][[j]]$smoothed_freq_growth <- predict(loess_freq_model)
#     
#   }
# }


smoothed.gr <- readRDS("Code/ProcessedData/Nov15_variants_gr.rds")


label.Korea <- c("South Korea-BA.2 (BA.2)", "South Korea-Q (B.1.1.7)")
label.UK <- c("United Kingdom-B.1", "United Kingdom-BA.5.2 (BA.5.2)")

pdf("Manuscripts/plots/S_G_smoothing.pdf", width = 10, height = 5) # Specify file name and size

par(mfrow = c(2, 2))

plot(
  x = data_list.days$days_30[["South Korea-BA.5.2 (BA.5.2)"]] %>%
    dplyr::select(gap) %>%
    unlist(),
  y = data_list.days$days_30[["South Korea-BA.5.2 (BA.5.2)"]] %>%
    dplyr::select(sharing) %>%
    unlist(),
     type = "p", pch = 3, cex = 0.6, xlab = "Time (Day)", ylab = "Share",
  main = "South Korea-BA.5.2 (BA.5.2)")

lines(
  x = data_list.days.smoothed$days_60[["South Korea-BA.5.2 (BA.5.2)"]] %>%
    dplyr::select(gap) %>%
    unlist(),
  y = data_list.days.smoothed$days_60[["South Korea-BA.5.2 (BA.5.2)"]] %>%
    dplyr::select(smoothed_sharing) %>%
    unlist(),
      type = "l", col = "blue", lwd = 2)

plot(
  x = smoothed.gr$days_30[[label.Korea[1]]] %>%
    dplyr::select(gap) %>%
    unlist(),
  y = smoothed.gr$days_30[[label.Korea[1]]] %>%
    dplyr::select(smoothed_sharing_growth) %>%
    unlist(),
  type = "p", pch = 6, cex = .5, col = "blue", xlab = "Time (Day)", ylab = "Share", 
  main = "South Korea-BA.2")


plot(
  x = data_list.days$days_30[[label.UK[2]]] %>%
    dplyr::select(gap) %>%
    unlist(),
  y = data_list.days$days_30[[label.UK[2]]] %>%
    dplyr::select(sharing) %>%
    unlist(),
  type = "p", pch = 3, cex = 0.6, xlab = "Time (Day)", ylab = "Share Growth Rate",
  main = "United Kingdom-BA.5.2")

lines(
  x = data_list.days.smoothed$days_60[[label.UK[2]]] %>%
    dplyr::select(gap) %>%
    unlist(),
  y = data_list.days.smoothed$days_60[[label.UK[2]]] %>%
    dplyr::select(smoothed_sharing) %>%
    unlist(),
  type = "l", col = "blue", lwd = 2)

plot(
  x = smoothed.gr$days_30[[label.UK[2]]] %>%
    dplyr::select(gap) %>%
    unlist(),
  y = smoothed.gr$days_30[[label.UK[2]]] %>%
    dplyr::select(smoothed_sharing_growth) %>%
    unlist(),
  type = "p", pch = 6, cex = .5, col = "blue", xlab = "Time (Day)", ylab = "Share Growth Rate", 
  main = "United Kingdom-BA.5.2")

dev.off()




# calculate the dominance? 

dominance <- data_list
for (i in names(data_list)) {
  
  # i = "USA-JN.1"
  
  data <- data_list[[i]]
  
  data <- data %>%
    mutate(
      group = with(rle(sharing != 0), rep(seq_along(values), lengths))
    ) %>%
    mutate(group = ifelse(sharing == 0, NA, group))  # Set group to NA for zero values
  
  # Filter out rows with zero sharing values and find the longest non-zero period
  longest_group <- data %>%
    dplyr::filter(!is.na(group)) %>%
    group_by(group) %>%
    summarise(length = n(), mean_sharing = mean(sharing)) %>%
    arrange(desc(length)) %>%
    slice(1)
  
  result <- matrix(nrow = 1, ncol = 2)
  colnames(result) <- c("country_label" ,"mean_share")
  
  result[1, 1] <- i
  result[1, 2] <- mean(longest_group$mean_sharing)
  
  dominance[[i]] <- result
  
}

dominance <- do.call(rbind, dominance) %>% as.data.frame()
dominance$mean_share <- as.numeric(dominance$mean_share)


start_duration <- matrix(nrow = dim(dominance)[1], ncol = 3) %>% as.data.frame()
colnames(start_duration) <- c("country_label", "start_date", "duration")
# get the duration and start date 
for (i in 1:length(names(data_list))) {
  
  # i = 1 # for test
  start_duration$country_label[i] <- names(data_list)[i]
  start_duration$start_date[i] <- data_list[[names(data_list)[i]]]$date[1] %>% as.character()
  start_duration$duration[i] <- nrow(data_list[[names(data_list)[i]]])
  
}


Nov15_variants_dominance.url <- "Code/ProcessedData/Nov15_variants_dominance.rds"
# saveRDS(dominance, Nov15_variants_dominance.url)


Nov15_variants_60days_gr.url <- "Code/ProcessedData/Nov15_variants_60days_gr.rds"
# saveRDS(smoothed.gr$days_60, Nov15_variants_60days_gr.url)

Nov25_variants_duration.url <- "Code/ProcessedData/Nov25_variants_duration.rds"
# saveRDS(start_duration, Nov15_variants_duration.url)



