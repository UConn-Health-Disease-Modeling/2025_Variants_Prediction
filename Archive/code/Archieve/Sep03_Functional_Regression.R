# load packages 
library(magrittr)
library(ggplot2)
library(fda)
source('Code/functions_franky_0725.R')

# load the data
url1 <- "Code/ProcessedData/list_30_.rds"
url2 <- "Code/ProcessedData/list_60_.rds"
url3 <- "Code/ProcessedData/list_dom_.rds"

# list.files()
list_30_ <- readRDS(url1)
list_60_ <- readRDS(url2)
list_dom <- readRDS(url3)

Korea_30 <- list_30_$Korea
Korea_60 <- list_60_$Korea
Korea_response <- list_dom$Korea

UK_30 <- list_30_$UK
UK_60 <- list_60_$UK
UK_response <- list_dom$UK


# convert to dataframes 
Korea_abs_raw_ls_ <- list()
Korea_share_raw_ls_ <- list()

Korea_abs_fit_ls_ <- list()
Korea_abs_gr_ls_ <- list()

Korea_share_fit_ls_ <- list()
Korea_share_gr_ls_ <- list()

for (name in names(Korea_30)) {
  
  # name <- "AY (B.1.617.2)"
  
  select_step <- seq(from = 1, to = 301, by = 10)
  
  abs_raw_col <- c(name, Korea_30[[name]]$raw$numerator %>% as.vector())
  share_raw_col <- c(name, Korea_30[[name]]$raw$sharing %>% as.vector())
  
  abs_fit_col <- c(name, (Korea_30[[name]]$fit$numerator_fit %>% as.vector())[select_step])
  abs_gr_col <- c(name, (Korea_30[[name]]$fit$numerator_gr_fit %>% as.vector())[select_step])
  
  share_fit_col <- c(name, (Korea_30[[name]]$fit$sharing_fit %>% as.vector())[select_step])
  share_gr_col <- c(name, (Korea_30[[name]]$fit$sharing_gr_fit %>% as.vector())[select_step])
  
  
  Korea_abs_raw_ls_[[name]] <- abs_raw_col
  Korea_share_raw_ls_[[name]] <- share_raw_col
  Korea_abs_fit_ls_[[name]] <- abs_fit_col
  Korea_abs_gr_ls_[[name]] <- abs_gr_col
  Korea_share_fit_ls_[[name]] <- share_fit_col
  Korea_share_gr_ls_[[name]] <- share_gr_col
  
}

Korea_abs_raw <- do.call(rbind, Korea_abs_raw_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
Korea_share_raw <- do.call(rbind, Korea_share_raw_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
Korea_abs_fit <- do.call(rbind, Korea_abs_fit_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
Korea_abs_gr <- do.call(rbind, Korea_abs_gr_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
Korea_share_fit <- do.call(rbind, Korea_share_fit_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
Korea_share_gr <- do.call(rbind, Korea_share_gr_ls_) %>% as.data.frame(stringsAsFactors = FALSE)


col_names <- c("classified_label", paste0("Day", 0:30))
colnames(Korea_abs_raw) <- col_names
colnames(Korea_share_raw) <- col_names
colnames(Korea_abs_fit) <- col_names
colnames(Korea_abs_gr) <- col_names
colnames(Korea_share_fit) <- col_names
colnames(Korea_share_gr) <- col_names

Korea_abs_raw <- Korea_abs_raw %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
Korea_share_raw <- Korea_share_raw %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
Korea_abs_fit <- Korea_abs_fit %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
Korea_abs_gr <- Korea_abs_gr %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
Korea_share_fit <- Korea_share_fit %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
Korea_share_gr <- Korea_share_gr %>% dplyr::mutate(across(starts_with("Day"), as.numeric))

rownames(Korea_abs_raw) <- NULL
rownames(Korea_share_raw) <- NULL
rownames(Korea_abs_fit) <- NULL
rownames(Korea_abs_gr) <- NULL
rownames(Korea_share_fit) <- NULL
rownames(Korea_share_gr) <- NULL

Korea_30_ls <- list(
  abs_row = Korea_abs_raw, 
  share_raw = Korea_share_raw, 
  abs_fit = Korea_abs_fit, 
  abs_gr = Korea_abs_gr, 
  share_fit = Korea_share_fit, 
  share_gr = Korea_share_gr
)


# UK 30 days
UK_abs_raw_ls_ <- list()
UK_share_raw_ls_ <- list()

UK_abs_fit_ls_ <- list()
UK_abs_gr_ls_ <- list()

UK_share_fit_ls_ <- list()
UK_share_gr_ls_ <- list()

for (name in names(UK_30)) {
  
  # name <- "AY (B.1.617.2)"
  
  select_step <- seq(from = 1, to = 301, by = 10)
  
  abs_raw_col <- c(name, UK_30[[name]]$raw$numerator %>% as.vector())
  share_raw_col <- c(name, UK_30[[name]]$raw$sharing %>% as.vector())
  
  abs_fit_col <- c(name, (UK_30[[name]]$fit$numerator_fit %>% as.vector())[select_step])
  abs_gr_col <- c(name, (UK_30[[name]]$fit$numerator_gr_fit %>% as.vector())[select_step])
  
  share_fit_col <- c(name, (UK_30[[name]]$fit$sharing_fit %>% as.vector())[select_step])
  share_gr_col <- c(name, (UK_30[[name]]$fit$sharing_gr_fit %>% as.vector())[select_step])
  
  
  UK_abs_raw_ls_[[name]] <- abs_raw_col
  UK_share_raw_ls_[[name]] <- share_raw_col
  UK_abs_fit_ls_[[name]] <- abs_fit_col
  UK_abs_gr_ls_[[name]] <- abs_gr_col
  UK_share_fit_ls_[[name]] <- share_fit_col
  UK_share_gr_ls_[[name]] <- share_gr_col
  
}

UK_abs_raw <- do.call(rbind, UK_abs_raw_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
UK_share_raw <- do.call(rbind, UK_share_raw_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
UK_abs_fit <- do.call(rbind, UK_abs_fit_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
UK_abs_gr <- do.call(rbind, UK_abs_gr_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
UK_share_fit <- do.call(rbind, UK_share_fit_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
UK_share_gr <- do.call(rbind, UK_share_gr_ls_) %>% as.data.frame(stringsAsFactors = FALSE)


col_names <- c("classified_label", paste0("Day", 0:30))
colnames(UK_abs_raw) <- col_names
colnames(UK_share_raw) <- col_names
colnames(UK_abs_fit) <- col_names
colnames(UK_abs_gr) <- col_names
colnames(UK_share_fit) <- col_names
colnames(UK_share_gr) <- col_names

UK_abs_raw <- UK_abs_raw %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
UK_share_raw <- UK_share_raw %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
UK_abs_fit <- UK_abs_fit %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
UK_abs_gr <- UK_abs_gr %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
UK_share_fit <- UK_share_fit %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
UK_share_gr <- UK_share_gr %>% dplyr::mutate(across(starts_with("Day"), as.numeric))

rownames(UK_abs_raw) <- NULL
rownames(UK_share_raw) <- NULL
rownames(UK_abs_fit) <- NULL
rownames(UK_abs_gr) <- NULL
rownames(UK_share_fit) <- NULL
rownames(UK_share_gr) <- NULL

UK_30_ls <- list(
  abs_row = UK_abs_raw, 
  share_raw = UK_share_raw, 
  abs_fit = UK_abs_fit, 
  abs_gr = UK_abs_gr, 
  share_fit = UK_share_fit, 
  share_gr = UK_share_gr
)





# for (name in names(Korea_60)) {
#   
#   # name <- "AY (B.1.617.2)"
#   
#   cat("raw sequence length: ", length(Korea_60[[name]][["raw"]]$numerator), "\n")
# }

korea_names <- c()
i = 1
for (name in names(Korea_60)) {
  if(length(Korea_60[[name]][["raw"]]$numerator) == 61){
    korea_names[i] <- name
    i = i + 1
  }
}



# convert to dataframes 
Korea_abs_raw_ls_ <- list()
Korea_share_raw_ls_ <- list()

Korea_abs_fit_ls_ <- list()
Korea_abs_gr_ls_ <- list()

Korea_share_fit_ls_ <- list()
Korea_share_gr_ls_ <- list()

for (name in korea_names) {
  
  # name <- "AY (B.1.617.2)"
  
  select_step <- seq(from = 1, to = 601, by = 10)
  
  abs_raw_col <- c(name, Korea_60[[name]]$raw$numerator %>% as.vector())
  share_raw_col <- c(name, Korea_60[[name]]$raw$sharing %>% as.vector())
  
  abs_fit_col <- c(name, (Korea_60[[name]]$fit$numerator_fit %>% as.vector())[select_step])
  abs_gr_col <- c(name, (Korea_60[[name]]$fit$numerator_gr_fit %>% as.vector())[select_step])
  
  share_fit_col <- c(name, (Korea_60[[name]]$fit$sharing_fit %>% as.vector())[select_step])
  share_gr_col <- c(name, (Korea_60[[name]]$fit$sharing_gr_fit %>% as.vector())[select_step])
  
  
  Korea_abs_raw_ls_[[name]] <- abs_raw_col
  Korea_share_raw_ls_[[name]] <- share_raw_col
  Korea_abs_fit_ls_[[name]] <- abs_fit_col
  Korea_abs_gr_ls_[[name]] <- abs_gr_col
  Korea_share_fit_ls_[[name]] <- share_fit_col
  Korea_share_gr_ls_[[name]] <- share_gr_col
  
}

Korea_abs_raw <- do.call(rbind, Korea_abs_raw_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
Korea_share_raw <- do.call(rbind, Korea_share_raw_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
Korea_abs_fit <- do.call(rbind, Korea_abs_fit_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
Korea_abs_gr <- do.call(rbind, Korea_abs_gr_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
Korea_share_fit <- do.call(rbind, Korea_share_fit_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
Korea_share_gr <- do.call(rbind, Korea_share_gr_ls_) %>% as.data.frame(stringsAsFactors = FALSE)


col_names <- c("classified_label", paste0("Day", 0:60))
colnames(Korea_abs_raw) <- col_names
colnames(Korea_share_raw) <- col_names
colnames(Korea_abs_fit) <- col_names
colnames(Korea_abs_gr) <- col_names
colnames(Korea_share_fit) <- col_names
colnames(Korea_share_gr) <- col_names

Korea_abs_raw <- Korea_abs_raw %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
Korea_share_raw <- Korea_share_raw %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
Korea_abs_fit <- Korea_abs_fit %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
Korea_abs_gr <- Korea_abs_gr %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
Korea_share_fit <- Korea_share_fit %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
Korea_share_gr <- Korea_share_gr %>% dplyr::mutate(across(starts_with("Day"), as.numeric))

rownames(Korea_abs_raw) <- NULL
rownames(Korea_share_raw) <- NULL
rownames(Korea_abs_fit) <- NULL
rownames(Korea_abs_gr) <- NULL
rownames(Korea_share_fit) <- NULL
rownames(Korea_share_gr) <- NULL

Korea_60_ls <- list(
  abs_row = Korea_abs_raw, 
  share_raw = Korea_share_raw, 
  abs_fit = Korea_abs_fit, 
  abs_gr = Korea_abs_gr, 
  share_fit = Korea_share_fit, 
  share_gr = Korea_share_gr
)


# UK 60 days
UK_abs_raw_ls_ <- list()
UK_share_raw_ls_ <- list()

UK_abs_fit_ls_ <- list()
UK_abs_gr_ls_ <- list()

UK_share_fit_ls_ <- list()
UK_share_gr_ls_ <- list()

for (name in names(UK_60)) {
  
  # name <- "AY (B.1.617.2)"
  
  select_step <- seq(from = 1, to = 601, by = 10)
  
  abs_raw_col <- c(name, UK_60[[name]]$raw$numerator %>% as.vector())
  share_raw_col <- c(name, UK_60[[name]]$raw$sharing %>% as.vector())
  
  abs_fit_col <- c(name, (UK_60[[name]]$fit$numerator_fit %>% as.vector())[select_step])
  abs_gr_col <- c(name, (UK_60[[name]]$fit$numerator_gr_fit %>% as.vector())[select_step])
  
  share_fit_col <- c(name, (UK_60[[name]]$fit$sharing_fit %>% as.vector())[select_step])
  share_gr_col <- c(name, (UK_60[[name]]$fit$sharing_gr_fit %>% as.vector())[select_step])
  
  
  UK_abs_raw_ls_[[name]] <- abs_raw_col
  UK_share_raw_ls_[[name]] <- share_raw_col
  UK_abs_fit_ls_[[name]] <- abs_fit_col
  UK_abs_gr_ls_[[name]] <- abs_gr_col
  UK_share_fit_ls_[[name]] <- share_fit_col
  UK_share_gr_ls_[[name]] <- share_gr_col
  
}

UK_abs_raw <- do.call(rbind, UK_abs_raw_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
UK_share_raw <- do.call(rbind, UK_share_raw_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
UK_abs_fit <- do.call(rbind, UK_abs_fit_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
UK_abs_gr <- do.call(rbind, UK_abs_gr_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
UK_share_fit <- do.call(rbind, UK_share_fit_ls_) %>% as.data.frame(stringsAsFactors = FALSE)
UK_share_gr <- do.call(rbind, UK_share_gr_ls_) %>% as.data.frame(stringsAsFactors = FALSE)


col_names <- c("classified_label", paste0("Day", 0:60))
colnames(UK_abs_raw) <- col_names
colnames(UK_share_raw) <- col_names
colnames(UK_abs_fit) <- col_names
colnames(UK_abs_gr) <- col_names
colnames(UK_share_fit) <- col_names
colnames(UK_share_gr) <- col_names

UK_abs_raw <- UK_abs_raw %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
UK_share_raw <- UK_share_raw %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
UK_abs_fit <- UK_abs_fit %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
UK_abs_gr <- UK_abs_gr %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
UK_share_fit <- UK_share_fit %>% dplyr::mutate(across(starts_with("Day"), as.numeric))
UK_share_gr <- UK_share_gr %>% dplyr::mutate(across(starts_with("Day"), as.numeric))

rownames(UK_abs_raw) <- NULL
rownames(UK_share_raw) <- NULL
rownames(UK_abs_fit) <- NULL
rownames(UK_abs_gr) <- NULL
rownames(UK_share_fit) <- NULL
rownames(UK_share_gr) <- NULL

UK_60_ls <- list(
  abs_row = UK_abs_raw, 
  share_raw = UK_share_raw, 
  abs_fit = UK_abs_fit, 
  abs_gr = UK_abs_gr, 
  share_fit = UK_share_fit, 
  share_gr = UK_share_gr
)

rm(list = setdiff(ls(), c("Korea_response", "UK_response", 
                          "Korea_30_ls", "UK_30_ls", "Korea_60_ls", "UK_60_ls", "paa_transform")))


################################################################################

# ls <- Korea_30_ls; response <- Korea_response # for test 
ls <- UK_30_ls; response <- UK_response # for test 

names_ <- names(ls)

output_df_ <- data.frame(
  Type = names_, 
  Rsq = numeric(length(names_))
)

for ( i in 1:length(names_)) {
  
  i = 4 # for test
  
  df_ <- ls[[names_[i]]] |> 
    dplyr::left_join(
      response, by = "classified_label"
    ) 
  
  # get the matrix
  input_ <- df_[2:32] |> 
    as.matrix() |> 
    t() 
  
  # get the output 
  y_ <- df_$dominance
  
  
  # generate the basis 
  smallbasis <- create.bspline.basis(c(0, 31), nbasis = 11)
  time_points <- seq(0, 31, length.out = 100)
  day.5 <- seq(0.5, 30.5, 1)
  
  tempfd <- smooth.basis(day.5, input_, smallbasis)$fd
  
  # fit the regression model 
  fRegress_model <- fRegress(y_ ~ tempfd, method = "fRegress")
  
  # coefficients_list <- fRegress_model$betaestlist
  # coef1 <- coefficients_list$const$fd$coefs
  # coef2 <- coefficients_list$tempfd$fd$coefs
  
  # get the fitted values 
  fit_ <- fRegress_model$yhatfdobj
  
  # calculate residuals
  residuals <- y_ - fit_
  
  # calculate tss and rss
  tss <- sum((y_ - mean(y_))^2)
  rss <- sum(residuals^2)
  
  # calculate R^2 
  r_squared <- 1 - (rss / tss)
  
  output_df_$Rsq[i] <- r_squared
  
}




















