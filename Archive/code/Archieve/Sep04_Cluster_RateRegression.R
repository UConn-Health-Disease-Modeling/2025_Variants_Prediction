# load the data 
url <- "Data/Sep04_combined_list.rds"
ls <- readRDS(url)

# source my function
source('Code/functions_franky_0725.R')

# names(ls)
# [1] "Korea_30_ls"    "Korea_60_ls"    "Korea_response" "UK_30_ls"       "UK_60_ls"       "UK_response" 

# it make sense to cluster by raw share 
Korea_share <- ls$Korea_30_ls$share_raw
UK_share <- ls$UK_30_ls$share_raw

# convert to list 
df2ls <- function(df){
  
  # df <- Korea_share # for test
  
  ls <- setNames(
    lapply(split(df[, -1], seq(nrow(df))), unlist),  
    df[, 1]                       
  )
  return(ls)
  
}

Korea_share_ls <- df2ls(Korea_share)
UK_share_ls <- df2ls(UK_share)


# install.packages("tsfeatures")
library(tsfeatures)

get_cluster <- function(ls, num_center = 3){
  
  # ls = Korea_share_ls # for test
  
    features <- tsfeatures(ls, features = c("mean", "var", "stl_features"))
  return(
    data.frame(
      classified_label = names(ls), 
      clusters = kmeans(features, centers = 3)$cluster
    )
  )
  
}

Korea_cluster <- get_cluster(Korea_share_ls)
UK_cluster <- get_cluster(UK_share_ls)

# more interested in UK performance since it is not good with original data

# 30 days 
len <- 3 * length(names(ls$UK_30_ls))

UK_result <- data.frame(
  country = rep("United Kingdom", len), 
  input = character(len), 
  cluster = rep(c(1, 2, 3), len/3), 
  Rsq = numeric(len), 
  pvalue = numeric(len)
)



k = 0

for (i in 1:length(names(ls$UK_30_ls))) {
  
  colnames_ <- c('response', '1-5', '6-10', '11-15', '16-20', '21-25', '26-30')
  
  # i = 1 # for test
  
  df_ <- ls$UK_30_ls[[names(ls$UK_30_ls)[i]]] |>
    left_join(ls$UK_response, by = "classified_label") |>
    left_join(UK_cluster, by = "classified_label")
  
  split_ls <- split(df_, df_$clusters)
  
  for (j in 1:length(names(split_ls))) {
    
    # j = 1 # for test 
    
    k = k + 1
    
    df_s_ <- split_ls[[names(split_ls)[j]]]
    
    X_ <- t(as.matrix(df_s_[, paste0("Day", 0:30)])) |>
      apply(2, paa_transform, num_segments = 6) |>
      t()
    
    Y_ <- df_s_$dominance
    
    data_input <- cbind(Y_, X_)
    colnames(data_input) <- colnames_
    
    model_ <- lm(response ~ ., data = data_input |> as.data.frame())
    
    UK_result$Rsq[k] <- summary(model_)$r.squared
    
    p_value <- summary(model_)$fstatistic
    p_value <- pf(p_value[1], p_value[2], p_value[3], lower.tail = FALSE)
    
    UK_result$pvalue[k] <- p_value
    UK_result$input[k] <- names(ls$UK_30_ls)[i]
    
  }
  
}

# visualize the clustering 
library(tidyr)
library(ggplot2)

plot_df <- ls$UK_30_ls$share_raw |>
  left_join(UK_cluster, by = "classified_label")


plot_long <- plot_df |>
  pivot_longer(cols = -c(classified_label, clusters), 
               names_to = "time",      
               values_to = "value")    

plot_long$time <- as.numeric(gsub("Day", "", colnames(plot_df)[2:32])) |> 
  rep(dim(plot_df)[1])


ggplot(plot_filtered, aes(x = time, y = value, group = classified_label, color = factor(clusters))) +
  geom_line() +
  labs(x = "Time", y = "Value", color = "ID") +
  facet_wrap(~ clusters) +
  theme_bw()


UK_result |> dplyr::filter(cluster == 2)
UK_result |> dplyr::filter(cluster == 1)
UK_result |> dplyr::filter(cluster == 3)
