library(dplyr)
# install.packages("glasso")
# install.packages("genlasso")
library(glasso)
library(genlasso)

mean_dominance_list.url <- "Code/ProcessedData/UKHSA_mean_dominance.rds"
save_url_rate.list.r <- "Code/ProcessedData/UKHSA_rate.list.r.rds"

mean_dominance_list <- readRDS(mean_dominance_list.url)
rate.list.r         <- readRDS(save_url_rate.list.r)

# names(mean_dominance_list)
for (country in names(mean_dominance_list)) {
  
  mean_dominance.country <- mean_dominance_list[[country]]
  
  for (i in names(rate.list.r)) {
    
    for (j in names(rate.list.r[[i]])) {
      
      rate.list.r[[i]][[j]][[country]] <- rate.list.r[[i]][[j]][[country]] |> 
        left_join(mean_dominance_list[[country]], by = "classified_label") 
      
    }
    
  }
  
}

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

# set the input
choice <- "sharing_fit.r"
names(data.60)


input_series <- data.60[[choice]] |> dplyr::select(starts_with("fit"))
output_dominance <- data.60[[choice]]$dominance
output_mean_dominance <- data.60[[choice]]$mean_dominance

# standardized_series <- t(apply(input_series, 1, function(x) {
#   (x - mean(x)) / sd(x)
# }))
# 
# standardized_series <- as.data.frame(standardized_series)


################################################################################
## convert to a occupational time curve 

seq.span <- input_series %>% as.matrix %>% as.vector()

bounds <- quantile(seq.span, probs = c(0.18, 0.98))
grid.span <- seq(from = bounds[1], to = bounds[2], length.out = 1000)

# output_dominance
input_otc <- input_series
input_otc[ , ] <- 0

for (i in 1:dim(input_series)[1]) {
  
  series <- input_series[i, ] %>% unlist %>% unname()
  # plot(series, ylim = c(bounds[1], bounds[2]))
  series.len <- length(series)
  
  series.cdf <- c()
  
  for (j in 1:length(grid.span)) {
    
    series.cdf[j] <- sum(series > grid.span[j])/series.len
    
  }
  
  # plot(x = grid.span, y = series.cdf)
  input_otc[i, ] <- series.cdf
  
}

# convert to a matrix 
input_otc <- input_otc %>% as.matrix() %>% unname()


# let the J to be 60 
matrix.A <- matrix(0, nrow = dim(input_otc)[1], ncol = 60) 

for (i in 1:60) {
  matrix.A[, i] <- rowSums(input_otc[, ((i - 1) * 10 + 1):(i * 10)])
}



fused_lasso_fit <- fusedlasso1d(y = train_y, X = train_X)

r_squared_values <- numeric(length(fused_lasso_fit$lambda))

for (i in 1:length(fused_lasso_fit$lambda)) {
  beta_hat <- coef(fused_lasso_fit)$beta[, i]  # Coefficients for current lambda
  predictions <- validation_X %*% beta_hat  # Predictions on the validation set
  
  SS_res <- sum((validation_y - predictions)^2)
  
  SS_tot <- sum((validation_y - mean(validation_y))^2)
  
  r_squared_values[i] <- 1 - (SS_res / SS_tot)
}










