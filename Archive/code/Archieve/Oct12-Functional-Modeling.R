library(dplyr)
library(fda)



mean_dominance_list.url <- "Code/ProcessedData/UKHSA_mean_dominance.rds"
# saveRDS(mean_dominance_list, mean_dominance_list.url)
save_url_rate.list.r <- "Code/ProcessedData/UKHSA_rate.list.r.rds"


mean_dominance_list <- readRDS(mean_dominance_list.url)
rate.list.r         <- readRDS(save_url_rate.list.r)

for (country in names(mean_dominance_list)) {
  
  # country = "USA"
  mean_dominance.country <- mean_dominance_list[[country]]
  
  for (i in names(rate.list.r)) {
    
    # i = "list.30"
    
    for (j in names(rate.list.r[[i]])) {
      
      # j = "sharing_fit"
      
      rate.list.r[[i]][[j]][[country]] <- rate.list.r[[i]][[j]][[country]] |> 
        left_join(mean_dominance_list[[country]], by = "classified_label") 
      
    }
    
  }
  
}

################################################################################
######################### 30 DAYS MODELING #####################################
################################################################################
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
output_mean_dominance <- data.30[[choice]]$mean_dominance

output_dominance_std <- (output_dominance - mean(output_dominance))/sd(output_dominance)
output_mean_dominance_std <- (output_mean_dominance - mean(output_mean_dominance))/sd(output_mean_dominance)

input_matrix <- input_series |> t() |> as.matrix()

# define the range 
rangeval <- c(0.1, 30)
simple_basis <- create.fourier.basis(rangeval,25)
day.index <- seq(0.1, by = 0.1, length.out = 300)

# fit the model 
Smooth <-  smooth.basis(day.index, input_matrix, simple_basis)
Smoothfd <- Smooth$fd

# plot(Smoothfd)
dominance.mdl <- fRegress(output_mean_dominance_std ~ Smoothfd)
plot(dominance.mdl$betaestlist[[2]]$fd)

# make prediction
dominance.pred = dominance.mdl$yhatfdobj
dominance.pred[dominance.pred<0] = 0

dominance.res = output_dominance_std - dominance.pred
SSE1  = sum(dominance.res^2)
SSTO    = sum((output_dominance_std - mean(output_dominance_std))^2)

Rsq <- 1 - SSE1/SSTO
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
  
  index_ = (data.30$sharing_fit.r$country == i)
  res_ <- output_dominance_std[index_] - dominance.pred[index_]
  sse <- sum(res_^2)
  ssto = sum((output_dominance_std[index_] - mean(output_dominance_std[index_]))^2)
  
  Rsq.list[[i]] = 1 - sse/ssto
  
  n <- length(output_dominance_std[index_]) # number of data points
  k <- 1  
  
  t_stat <- sqrt(Rsq.list[[i]]) * sqrt((n - k - 1) / (1 - Rsq.list[[i]]))
  
  p_value <- 2 * pt(-abs(t_stat), df = n - k - 1)
  
  p_value.list[[i]] = p_value
}

Rsq_p_values.df.30 <- data.frame(
  country = names(Rsq.list),
  Rsq_value.30 = unlist(Rsq.list),
  p_value.30 = unlist(p_value.list)
) |> arrange(desc(Rsq_value.30))


# country_order <- Rsq_p_values.df.30$country
################################################################################
######################### END 30 DAYS MODELING #################################
################################################################################










################################################################################
######################### 45 DAYS MODELING #####################################
################################################################################
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
output_mean_dominance <- data.45[[choice]]$mean_dominance

output_dominance_std <- (output_dominance - mean(output_dominance))/sd(output_dominance)
output_mean_dominance_std <- (output_mean_dominance - mean(output_mean_dominance))/sd(output_mean_dominance)

input_matrix <- input_series |> t() |> as.matrix()

rangeval <- c(0.1, 45)
simple_basis <- create.fourier.basis(rangeval,21)

day.index <- seq(0.1, by = 0.1, length.out = 450)
Smooth <-  smooth.basis(day.index, input_matrix, simple_basis)
Smoothfd <- Smooth$fd

plot(Smoothfd)
dominance.mdl <- fRegress(output_mean_dominance_std ~ Smoothfd)
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
################################################################################
######################### END 45 DAYS MODELING #################################
################################################################################










################################################################################
######################### 60 DAYS MODELING #####################################
################################################################################
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
output_mean_dominance <- data.60[[choice]]$mean_dominance

output_dominance_std <- (output_dominance - mean(output_dominance))/sd(output_dominance)
output_mean_dominance_std <- (output_mean_dominance - mean(output_mean_dominance))/sd(output_mean_dominance)


rangeval <- c(0.1, 60)
simple_basis <- create.fourier.basis(rangeval,21)

day.index <- seq(0.1, by = 0.1, length.out = 600)
Smooth <-  smooth.basis(day.index, input_matrix, simple_basis)
Smoothfd <- Smooth$fd

plot(Smoothfd)
dominance.mdl <- fRegress(output_mean_dominance_std ~ Smoothfd)
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
################################################################################
######################### END 60 DAYS MODELING #################################
################################################################################












Rsq_p_values.df.comb <- Rsq_p_values.df.30 |> 
  left_join(Rsq_p_values.df.45, by = "country") |> 
  left_join(Rsq_p_values.df.60, by = "country")

df1 <- Rsq_p_values.df.comb |> dplyr::select(country, starts_with("Rsq"))
df2 <- Rsq_p_values.df.comb |> dplyr::select(country, starts_with("p_value"))


df2 <- df2 |> 
  mutate(across(starts_with("p_value"), 
                ~ case_when(
                  . < 0.001 ~ "<0.001 ***",
                  . >= 0.001 & . < 0.01 ~ "<0.01 **",
                  . >= 0.01 & . < 0.05 ~ "<0.05 *",
                  . >= 0.05 & . < 0.1 ~ "<0.1 .",
                  TRUE ~ ">0.1"
                )))

# df1 <- df1 |>
#   mutate(country = factor(country, levels = country_order)) |>
#   arrange(country)
# 
# df2 <- df2 |>
#   mutate(country = factor(country, levels = country_order)) |>
#   arrange(country)

write.csv(df1, "../../../Desktop/df1.csv", row.names = FALSE)
write.csv(df2, "../../../Desktop/df2.csv", row.names = FALSE)

