library(dplyr)
library(magrittr)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(viridis)
library(scales) 
library(fda)
library(fregion)
library(boot)
library(splines)

source('Background/functions_0706.R')
source('Code/functions_franky_0725.R')
source('Background/growth_rate_function_original.R')

url <- "UKHSA-UConn-variant-modelling/variant_modelling/data/summary_GISAID_20240918.csv"

data <- read.csv(url)

main_countries <- data |> 
  group_by(country) |> 
  summarise(num_total = sum(numerator)) |> 
  arrange(desc(num_total)) |> 
  head(30) |> 
  dplyr::select(country) |> unlist() |> as.vector()

data2 <- data |> 
  filter(country %in% main_countries)

length(unique(data2$lineage)) # 3,746


# classify the lineage within the countries
# get the alias
alias_list <- get_alias()

data3 <- data2 |> 
  rename(n = numerator) |> 
  dplyr::select(country, date, lineage, n)


# classified the variants
summarised <- data3 %>%
  dplyr::group_by(country, lineage) %>%
  dplyr::summarise(n = sum(n))

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
# 12499           8448          12865           9762          24821           4517           6831          17079 
# France        Germany          India        Ireland         Israel          Italy          Japan     Luxembourg 
# 17930          28132          10075           8534          10179          13593          16707           7491 
# Mexico    Netherlands         Norway           Peru         Poland         Russia       Slovenia    South Korea 
# 9618          14879           8657           5349           6935           8134           7426           9941 
# Spain         Sweden    Switzerland         Turkey United Kingdom            USA 
# 16245          14579          14514           3222          42522          72354 

# classified_data |> filter(country == "Australia")

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

data_list <- list()
dominance_list <- list()

for (i in unique(classified_data$country)) {
  
  # i = "South Korea" # for test
  
  country_data <- classified_data |> 
    filter(country == i) |> 
    arrange(classified_label, date) |> 
    group_by(classified_label) |> 
    filter(n() > 60) |> 
    mutate(
      sharing = dplyr::coalesce(share, 0),
      gap = as.numeric(date - dplyr::first(date))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-share)
  
  data_list[[i]] <- country_data
  dominance_list[[i]] <- country_data %>% 
    dplyr::group_by(
      classified_label
    ) %>% 
    dplyr::summarise(
      dominance = sum(sharing, na.rm = TRUE),
    ) %>% 
    dplyr::arrange(
      desc(dominance)
    )
  
}

data_list.30 <- list()
data_list.45 <- list()
data_list.60 <- list()

for (i in names(data_list)) {
  
  data_list.30[[i]] <- data_list[[i]] |> 
    filter(gap <= 30)
  
  data_list.45[[i]] <- data_list[[i]] |> 
    filter(gap <= 45)
  
  data_list.60[[i]] <- data_list[[i]] |> 
    filter(gap <= 60)
  
}



################################################################################
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
################################################################################



rate_list.30 <- list()
rate_list.45 <- list()
rate_list.60 <- list()

for (i in names(data_list.30)) {
  
  rate_list.30[[i]] <- generate_growth_rate_list(data_list.30[[i]])
  rate_list.45[[i]] <- generate_growth_rate_list(data_list.45[[i]])
  rate_list.60[[i]] <- generate_growth_rate_list(data_list.60[[i]])
  
}

# save the result 
rate.list <- list(
  # list.30.raw = data_list.30, 
  # list.45.raw = data_list.45, 
  # list.60.raw = data_list.60, 
  list.30 = rate_list.30, 
  list.45 = rate_list.45, 
  list.60 = rate_list.60
)

# names(rate.list[["list.30"]])

save_url <- "Code/ProcessedData/UKHSA_rate_list.rds"
# saveRDS(rate.list, save_url)



plot_variant_data <- function(rate_list, country_name) {
  variants_list <- rate_list[[country_name]]
  
  variant_colors <- setNames(sample(colors(), length(variants_list)), names(variants_list))
  
  par(mar = c(5, 4, 4, 8) + 0.1)
  
  plot(NULL, xlim = c(1, length(variants_list[[1]]$raw$sharing)), 
       ylim = range(unlist(lapply(variants_list, function(x) x$raw$sharing))),
       xlab = "Index", ylab = "Raw Value", main = paste("Raw Data for All", country_name, "Variants"))
  
  for (variant_name in names(variants_list)) {
    raw_data <- variants_list[[variant_name]]$raw$sharing
    lines(raw_data, type = "l", lwd = 2, col = variant_colors[variant_name])
  }
  
  legend("topright", inset = c(-0.3, 0), legend = names(variants_list), 
         col = variant_colors, lty = 1, cex = 0.6, xpd = TRUE)
  
  par(mar = c(5, 4, 4, 8) + 0.1)
  
  plot(NULL, xlim = c(1, length(variants_list[[1]]$fit$sharing_fit)), 
       ylim = range(unlist(lapply(variants_list, function(x) x$fit$sharing_fit))),
       xlab = "Index", ylab = "Fitted Value", main = paste("Fitted Data for All", country_name, "Variants"))
  
  for (variant_name in names(variants_list)) {
    raw_data <- variants_list[[variant_name]]$fit$sharing_fit
    lines(raw_data, type = "l", lwd = 2, col = variant_colors[variant_name])
  }
  
  legend("topright", inset = c(-0.3, 0), legend = names(variants_list), 
         col = variant_colors, lty = 1, cex = 0.6, xpd = TRUE)
}

# Example usage:
plot_variant_data(rate_list.60, "Australia")
plot_variant_data(rate_list.60, "USA")
plot_variant_data(rate_list.60, "United Kingdom")
plot_variant_data(rate_list.60, "South Korea")



################################################################################
# fit the model with rate_list.30 
rate.list <- readRDS(save_url)

# rate.list$list.30$Australia$AY.103$fit |> colnames()
# "Date_numeric"     "date"             "numerator_fit"    "numerator_gr_fit" "sharing_fit"      "sharing_gr_fit" 
# Generalized function to handle rate list processing
process_rate_list <- function(rate_list, list_name, num_cols, dominance_list) {
  rate_list_r <- list(
    numerator_fit = list(),
    numerator_gr_fit = list(),
    sharing_fit = list(),
    sharing_gr_fit = list()
  )
  
  # Helper function to create dataframes for each type of fit
  create_fit_df <- function(country_list, fit_type, num_cols) {
    fit_df <- matrix(NA, nrow = length(names(country_list)), ncol = (1 + num_cols)) |> as.data.frame()
    for (j in 1:length(names(country_list))) {
      fit_data <- country_list[[j]]$fit[[fit_type]] |> unlist() |> as.vector()
      fit_df[j, 1] <- names(country_list)[j]
      fit_df[j, 2:(num_cols + 1)] <- fit_data
    }
    colnames(fit_df) <- c("classified_label", paste0("fit", seq(0, num_cols / 10, by = 0.1)))
    return(fit_df)
  }
  
  # Process each country list in the rate list
  for (i in names(rate_list)) {
    country_list <- rate_list[[i]]
    
    # Create dataframes for each fit type
    rate_list_r$numerator_fit[[i]] <- create_fit_df(country_list, "numerator_fit", num_cols)
    rate_list_r$numerator_gr_fit[[i]] <- create_fit_df(country_list, "numerator_gr_fit", num_cols)
    rate_list_r$sharing_fit[[i]] <- create_fit_df(country_list, "sharing_fit", num_cols)
    rate_list_r$sharing_gr_fit[[i]] <- create_fit_df(country_list, "sharing_gr_fit", num_cols)
  }
  
  # Apply left_join and select operations
  apply_dominance_join <- function(fit_list, dominance_list) {
    for (i in names(fit_list)) {
      fit_list[[i]] <- fit_list[[i]] |>
        left_join(dominance_list[[i]], by = "classified_label") |> 
        dplyr::select(-fit0)
    }
    return(fit_list)
  }
  
  # Apply dominance join for each fit type
  rate_list_r$numerator_fit <- apply_dominance_join(rate_list_r$numerator_fit, dominance_list)
  rate_list_r$numerator_gr_fit <- apply_dominance_join(rate_list_r$numerator_gr_fit, dominance_list)
  rate_list_r$sharing_fit <- apply_dominance_join(rate_list_r$sharing_fit, dominance_list)
  rate_list_r$sharing_gr_fit <- apply_dominance_join(rate_list_r$sharing_gr_fit, dominance_list)
  
  return(rate_list_r)
}

# Process each rate list (list.30, list.45, list.60)
rate.list.r <- list(
  list.30 = process_rate_list(rate.list[["list.30"]], "list.30", 301, dominance_list),
  list.45 = process_rate_list(rate.list[["list.45"]], "list.45", 451, dominance_list),
  list.60 = process_rate_list(rate.list[["list.60"]], "list.60", 601, dominance_list)
)


save_url_rate.list.r <- "Code/ProcessedData/UKHSA_rate.list.r.rds"
# saveRDS(rate.list.r, save_url_rate.list.r)






################################################################################
rate.list.r <- readRDS(save_url_rate.list.r)

# 30 days model
data.30 <- rate.list.r$list.30

for (i in names(data.30)) {

  # i <- "numerator_fit"

  rate.df <- data.30[[i]]

  for (j in 1:length(names(rate.df))){

    # j = 1

    if(j == 1){
      comb.df <- rate.df[[names(rate.df)[j]]] |>
        mutate(country = names(rate.df)[j])
    }else{
      comb.df <- comb.df |> rbind(rate.df[[names(rate.df)[j]]] |> mutate(country = names(rate.df)[j]))
    }
  }
  data.30[[paste0(i, ".r")]] <- comb.df
}

choice <- names(data.30)[7]
input_series <- data.30[[choice]] |> dplyr::select(starts_with("fit"))
output_dominance <- data.30[[choice]]$dominance

output_dominance_std <- (output_dominance - mean(output_dominance))/sd(output_dominance)

input_matrix <- input_series |> t() |> as.matrix()

rangeval <- c(0.1, 30)
simple_basis <- create.fourier.basis(rangeval,25)

day.index <- seq(0.1, by = 0.1, length.out = 300)
Smooth <-  smooth.basis(day.index, input_matrix, simple_basis)
Smoothfd <- Smooth$fd

plot(Smoothfd)
dominance.mdl <- fRegress(output_dominance_std ~ Smoothfd)
plot(dominance.mdl$betaestlist[[2]]$fd)


dominance.pred = dominance.mdl$yhatfdobj
dominance.pred[dominance.pred<0] = 0

dominance.res = output_dominance_std - dominance.pred
SSE1  = sum(dominance.res^2)

# sum squared residuals for the null model y = alpha + \epsilon
SSTO    = sum((output_dominance_std - mean(output_dominance_std))^2)

Rsq <- 1 - SSE1/SSTO
Rsq

Rsq.list <- list()
for (i in unique(data.30$sharing_fit.r$country)) {
  
  # i = "USA"
  index_ = (data.30$sharing_fit.r$country == i)
  res_ <- output_dominance_std[index_] - dominance.pred[index_]
  sse <- sum(res_^2)
  ssto    = sum((output_dominance_std[index_] - mean(output_dominance_std[index_]))^2)
  
  Rsq.list[[i]] = 1 - sse/ssto
}

p_value.list <- list()

for (i in unique(data.30$sharing_fit.r$country)) {
  
  # i = "USA"
  index_ = (data.30$sharing_fit.r$country == i)
  res_ <- output_dominance_std[index_] - dominance.pred[index_]
  sse <- sum(res_^2)
  ssto = sum((output_dominance_std[index_] - mean(output_dominance_std[index_]))^2)
  
  # Calculate R-squared
  Rsq.list[[i]] = 1 - sse/ssto
  
  # Calculate degrees of freedom
  n <- length(output_dominance_std[index_]) # number of data points
  k <- 1  # number of predictors (since we are using smoothfd as a single predictor)
  
  # Calculate the t-statistic (using the assumption that residuals are normally distributed)
  t_stat <- sqrt(Rsq.list[[i]]) * sqrt((n - k - 1) / (1 - Rsq.list[[i]]))
  
  # Calculate the p-value using the t-distribution
  p_value <- 2 * pt(-abs(t_stat), df = n - k - 1)
  
  p_value.list[[i]] = p_value
}

# Create a dataframe to combine the Rsq and p-values for each country
Rsq_p_values.df.30 <- data.frame(
  country = names(Rsq.list),
  Rsq_value.30 = unlist(Rsq.list),
  p_value.30 = unlist(p_value.list)
) |> arrange(desc(Rsq_value.30))




# 45 days model
data.45 <- rate.list.r$list.45

for (i in names(data.45)) {
  
  # i <- "numerator_fit"
  
  rate.df <- data.45[[i]]
  
  for (j in 1:length(names(rate.df))){
    
    # j = 1
    
    if(j == 1){
      comb.df <- rate.df[[names(rate.df)[j]]] |>
        mutate(country = names(rate.df)[j])
    }else{
      comb.df <- comb.df |> rbind(rate.df[[names(rate.df)[j]]] |> mutate(country = names(rate.df)[j]))
    }
  }
  data.45[[paste0(i, ".r")]] <- comb.df
}


choice <- names(data.45)[7]
input_series <- data.45[[choice]] |> dplyr::select(starts_with("fit"))
output_dominance <- data.45[[choice]]$dominance

output_dominance_std <- (output_dominance - mean(output_dominance))/sd(output_dominance)

input_matrix <- input_series |> t() |> as.matrix()

rangeval <- c(0.1, 45)
simple_basis <- create.fourier.basis(rangeval,21)

day.index <- seq(0.1, by = 0.1, length.out = 450)
Smooth <-  smooth.basis(day.index, input_matrix, simple_basis)
Smoothfd <- Smooth$fd

plot(Smoothfd)
dominance.mdl <- fRegress(output_dominance_std ~ Smoothfd)
plot(dominance.mdl$betaestlist[[2]]$fd)


dominance.pred = dominance.mdl$yhatfdobj
dominance.pred[dominance.pred<0] = 0

dominance.res = output_dominance_std - dominance.pred
SSE1  = sum(dominance.res^2)

# sum squared residuals for the null model y = alpha + \epsilon
SSTO    = sum((output_dominance_std - mean(output_dominance_std))^2)

Rsq <- 1 - SSE1/SSTO
Rsq

Rsq.list <- list()
for (i in unique(data.45$sharing_fit.r$country)) {
  
  # i = "USA"
  index_ = (data.45$sharing_fit.r$country == i)
  res_ <- output_dominance_std[index_] - dominance.pred[index_]
  sse <- sum(res_^2)
  ssto    = sum((output_dominance_std[index_] - mean(output_dominance_std[index_]))^2)
  
  Rsq.list[[i]] = 1 - sse/ssto
}


p_value.list <- list()

for (i in unique(data.45$sharing_fit.r$country)) {
  
  # i = "USA"
  index_ = (data.45$sharing_fit.r$country == i)
  res_ <- output_dominance_std[index_] - dominance.pred[index_]
  sse <- sum(res_^2)
  ssto = sum((output_dominance_std[index_] - mean(output_dominance_std[index_]))^2)
  
  # Calculate R-squared
  Rsq.list[[i]] = 1 - sse/ssto
  
  # Calculate degrees of freedom
  n <- length(output_dominance_std[index_]) # number of data points
  k <- 1  # number of predictors (since we are using smoothfd as a single predictor)
  
  # Calculate the t-statistic (using the assumption that residuals are normally distributed)
  t_stat <- sqrt(Rsq.list[[i]]) * sqrt((n - k - 1) / (1 - Rsq.list[[i]]))
  
  # Calculate the p-value using the t-distribution
  p_value <- 2 * pt(-abs(t_stat), df = n - k - 1)
  
  p_value.list[[i]] = p_value
}

# Create a dataframe to combine the Rsq and p-values for each country
Rsq_p_values.df.45 <- data.frame(
  country = names(Rsq.list),
  Rsq_value.45 = unlist(Rsq.list),
  p_value.45 = unlist(p_value.list)
) |> arrange(desc(Rsq_value.45))



# 60 days model
data.60 <- rate.list.r$list.60

for (i in names(data.60)) {
  
  # i <- "numerator_fit"
  
  rate.df <- data.60[[i]]
  
  for (j in 1:length(names(rate.df))){
    
    # j = 1
    
    if(j == 1){
      comb.df <- rate.df[[names(rate.df)[j]]] |>
        mutate(country = names(rate.df)[j])
    }else{
      comb.df <- comb.df |> rbind(rate.df[[names(rate.df)[j]]] |> mutate(country = names(rate.df)[j]))
    }
  }
  data.60[[paste0(i, ".r")]] <- comb.df
}


choice <- names(data.60)[7]
input_series <- data.60[[choice]] |> dplyr::select(starts_with("fit"))
output_dominance <- data.60[[choice]]$dominance

output_dominance_std <- (output_dominance - mean(output_dominance))/sd(output_dominance)

input_matrix <- input_series |> t() |> as.matrix()

rangeval <- c(0.1, 60)
simple_basis <- create.fourier.basis(rangeval,21)

day.index <- seq(0.1, by = 0.1, length.out = 600)
Smooth <-  smooth.basis(day.index, input_matrix, simple_basis)
Smoothfd <- Smooth$fd

plot(Smoothfd)
dominance.mdl <- fRegress(output_dominance_std ~ Smoothfd)
plot(dominance.mdl$betaestlist[[2]]$fd)


dominance.pred = dominance.mdl$yhatfdobj
dominance.pred[dominance.pred<0] = 0

dominance.res = output_dominance_std - dominance.pred
SSE1  = sum(dominance.res^2)

# sum squared residuals for the null model y = alpha + \epsilon
SSTO    = sum((output_dominance_std - mean(output_dominance_std))^2)

Rsq <- 1 - SSE1/SSTO
Rsq

Rsq.list <- list()
for (i in unique(data.60$sharing_fit.r$country)) {
  
  # i = "USA"
  index_ = (data.60$sharing_fit.r$country == i)
  res_ <- output_dominance_std[index_] - dominance.pred[index_]
  sse <- sum(res_^2)
  ssto    = sum((output_dominance_std[index_] - mean(output_dominance_std[index_]))^2)
  
  Rsq.list[[i]] = 1 - sse/ssto
}

p_value.list <- list()

for (i in unique(data.60$sharing_fit.r$country)) {
  
  # i = "USA"
  index_ = (data.60$sharing_fit.r$country == i)
  res_ <- output_dominance_std[index_] - dominance.pred[index_]
  sse <- sum(res_^2)
  ssto = sum((output_dominance_std[index_] - mean(output_dominance_std[index_]))^2)
  
  # Calculate R-squared
  Rsq.list[[i]] = 1 - sse/ssto
  
  # Calculate degrees of freedom
  n <- length(output_dominance_std[index_]) # number of data points
  k <- 1  # number of predictors (since we are using smoothfd as a single predictor)
  
  # Calculate the t-statistic (using the assumption that residuals are normally distributed)
  t_stat <- sqrt(Rsq.list[[i]]) * sqrt((n - k - 1) / (1 - Rsq.list[[i]]))
  
  # Calculate the p-value using the t-distribution
  p_value <- 2 * pt(-abs(t_stat), df = n - k - 1)
  
  p_value.list[[i]] = p_value
}

# Create a dataframe to combine the Rsq and p-values for each country
Rsq_p_values.df.60 <- data.frame(
  country = names(Rsq.list),
  Rsq_value.60 = unlist(Rsq.list),
  p_value.60 = unlist(p_value.list)
) |> arrange(desc(Rsq_value.60))


Rsq_p_values.df.comb <- Rsq_p_values.df.30 |> 
  left_join(Rsq_p_values.df.45, by = "country") |> 
  left_join(Rsq_p_values.df.60, by = "country")

Rsq_p_values.df.comb |> dplyr::select(country, starts_with("Rsq"))
Rsq_p_values.df.comb |> dplyr::select(country, starts_with("p_value"))
