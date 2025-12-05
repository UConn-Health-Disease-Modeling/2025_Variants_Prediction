#' Code for calculating variant growth rates, both relative and absolute. To 
#' apply this to all variants of interest, it needs to be looped over different
#' variants. 
#' 
#' Relative growth rate - this is the growth rate relative to another variant.
#' That is, if we look at the proportion of cases that are variant A relative to
#' the number of cases that are either variant A or variant B. Here, variant B 
#' can be a single variant (baseline = "variant_name") or every other variant 
#' (baseline = "all-other"). The all-other scenario captures the rate at which
#' variant A is taking over the whole population. The variant_name scenario
#' captures the rate at which variant A is replacing just variant B.  
#' 
#' Absolute growth rate - this is the growth rate for the number of cases that 
#' are variant A. 
#' 
#' Options:
#' 
#' denominator - this controls the smoothness of the fit. Larger value results
#' in increase smoothness. 
#' 
#' dow - this adjusts for day of week effects, options are "none" - no day of 
#' week effect, "re" - random effects forday of week, or "fixed" - fixed effects 
#' for day of week. I suggest "none" or "re". 


library(lubridate)
library(dplyr)
library(mgcv)
library(ggplot2)

setwd("D://My Drive//6. Korea Modeling//Variants//Code//")
df <- read.csv("variants_data_test.csv") # Load example data
source("growth_rate_function_original.R")

## Which variant are you interested in?
variant_ = "V-23JUL-01 (Omicron EG.5.1)" # "variant_name"

## What would you like to compare the relative growth to?
baseline = "V-23APR-01 (Omicron XBB.1.16)" # "variant_name" or "all-other" or "none"

## Would you like to calculate relative growth or absolute growth?
type = "absolute-growth" # "relative-growth" or "absolute-growth"

## Create input data for variant of interest
df_variant1 <- df %>% filter(variant == variant_)
date_min = df_variant1 %>% filter(n_samples>0) 
date_min = min(date_min$date)
df_variant1 <- df %>% filter(variant == variant_,date >= date_min)
N = df_variant1$n_samples
date_ = df_variant1$date

## Create input data for comparison data
if (baseline == "all-other"){
  offset = df_variant1$weekly_seq
} else if (baseline == "none"){
  offset = NA
} else {
  df_variant2 <- df%>% 
    filter(variant == baseline,date >= date_min)
  offset = df_variant2$n_samples + N + 0.00001
}


## Run the appropriate type of model
if (type == "relative-growth"){
  output <- extract_growth_rates(N,date_,offset,dow = "none",denominator = 7)
}
if (type == "absolute-growth"){
  output <- extract_growth_rates(N,date_,dow = "none",denominator = 7)
}


## Plot growth rate
growth_rates <- output$fitted_rates
ggplot(growth_rates,aes(x=date)) +
  geom_line(aes(y = rate),color = "red") + 
  geom_ribbon(aes(ymin = lower_rate, ymax = upper_rate),fill = "red", alpha = 0.1) 

## Plot model fit
fitted_model <- output$model_fit %>% mutate(date = as.Date(date))
if (type == "relative-growth"){
  plotting_data <- df_variant1 %>%
    mutate(y_var = n_samples/offset)
}
if (type == "absolute-growth"){
  plotting_data <- df_variant1 %>%
    mutate(y_var = n_samples)
}
ggplot(fitted_model,aes(x=date)) +
  geom_line(aes(y = full_fit),color = "blue") + 
  geom_ribbon(aes(ymin = full_lower, ymax = full_upper),fill = "blue",alpha = 0.1) +
  geom_point(data = plotting_data,aes(x = as.Date(date), y= y_var))
