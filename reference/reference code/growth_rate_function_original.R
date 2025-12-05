
require(mgcv)
require(dplyr)

extract_growth_rates <- function(N,
                                 date,
                                 offset = NULL,
                                 denominator = 10,
                                 bs = 'cr',
                                 dow = c('none', 're', 'fixed'),
                                 fineness = 0.1){

  if(length(N) != length(date)){stop('N and date time series are unequal lengths')}

  date <- as_date(date)

  Date_numeric <- as.numeric(date - min(date))
  d <- as.factor(weekdays(date))
  this_k <-  floor(length(Date_numeric) / denominator)

  the_formula <- switch(dow,
                        none = as.formula(N ~ s(Date_numeric, k = this_k, bs = bs)),
                        re = as.formula(N ~ s(Date_numeric, k = this_k, bs = bs) + s(d, bs = 're')),
                        fixed = as.formula(N ~ s(Date_numeric, k = this_k, bs = bs) + d))


  if(is.null(offset)){

    this_model <- gam(formula = the_formula,
                      family = 'nb')


  }else{

    this_model <- gam(formula = the_formula,
                      offset = log(offset),
                      family = 'nb')

  }


  data.frame(Date_numeric = seq(min(Date_numeric), max(Date_numeric), fineness)) %>%
    mutate(date = Date_numeric + min(date),
           d = as.factor(weekdays(date))) -> new_data0

  full_predictions <- predict(object = this_model,
                              type = 'link',
                              se.fit = TRUE,
                              newdata = new_data0,
                              exclude = 's(d)')

  new_data0 %>%
    mutate(full_fit = exp(full_predictions$fit),
           full_upper = exp(full_predictions$fit + 2 * full_predictions$se.fit),
           full_lower = exp(full_predictions$fit - 2 * full_predictions$se.fit))  -> model_fit

  matrix0 <- predict(object = this_model,
                     newdata = new_data0,
                     type = 'lpmatrix',
                     exclude = 's(d)')

  new_data0 %>%
    mutate(Date_numeric = Date_numeric + 1e-7) -> new_data1

  matrix1 <- predict(object = this_model,
                     newdata = new_data1,
                     type = 'lpmatrix',
                     exclude = 's(d)')

  X <- (matrix1 - matrix0) / 1e-7

  new_data0 %>%
    mutate(rate = as.vector(X %*% coef(this_model)),
           upper_rate = as.vector(rate + 2 * rowSums(X %*% this_model$Vp * X) ^ 0.5),
           lower_rate = as.vector(rate - 2 * rowSums(X %*% this_model$Vp * X) ^ 0.5)) -> fitted_rates

  return(list(fitted_rates = fitted_rates,
              model_fit = model_fit,
              the_model = this_model))

}
