# source the functions (functions_0706.R)
library(dplyr)
source(glue::glue('/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants_Analysis/Background/functions_0706.R'))
source('Code/functions_franky_0725.R')

################################################################################
# DON'T RUN
# list.files()
url <- "gisaid_variants_statistics_2024_08_05_1942/gisaid_variants_statistics.json"
data <- jsonlite::fromJSON(url)$stats

names(data)

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
    
    to_append <- country_of_date$submissions_per_clade
    
    if(length(to_append) > 0){
      
      to_append$date <- date
      to_append$country <- country
      
      data_df <- rbind(data_df, to_append)
      
    } 
    
  }
  
}

colnames(data_df) <- c("count", "clade", "date", "country")

data_df <- data_df |> 
  select(
    country, date, clade, count
  )

# save the dataset
save_url <- "Data/gisaid_variants_statistics.rds"
saveRDS(data_df, save_url)
################################################################################

data_df |> filter(country == "South Korea")

length(unique(data_df$clade))











