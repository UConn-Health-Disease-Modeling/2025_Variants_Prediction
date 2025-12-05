
# install packages
source('Collab/functions/setup.R')

source('Collab/functions/functions.R')

# summary pull from GISAID: https://gisaid.org/hcov19-variants/
dataset <- readr::read_csv('Collab/data/summary_gisaid.csv')

# Each row is a country , day, lineage, and count of total tests sequenced (denominator) and those for that lineage (numerator)
head(dataset)

# Test run - compare BQ.1 lineages for SK and UK 
lineage <- 'XBB.1.5'

# No alias needed here

# alias_list <- get_alias()

# this_alias <- dplyr::filter(alias_list, alias == 'XBB') |>
#   dplyr::pull(alias_lineage); this_alias

# append alias to data set
with_alias <- dataset |>
  dplyr::mutate(
    alias = gsub(
      pattern = "[^a-zA-Z]",
      replacement = "",
      x = lineage 
    ),
    second_part = gsub(pattern = "[a-zA-Z]",
                       replacement = "",
                       x = lineage)
  ) |>
  dplyr::left_join(alias_list) |>
  dplyr::mutate(
    unaliased_lineage  =
      dplyr::if_else(!is.na(alias_lineage),
                     paste0(alias_lineage, second_part),
                     lineage),
  )

summarised <- with_alias |>
  dplyr::mutate(
    is_target_lineage = stringr::str_starts(unaliased_lineage, 'XBB.1.5.') | unaliased_lineage == 'XBB.1.5'
  ) |>
  dplyr::summarise(
    numerator = sum(n * is_target_lineage), 
    denominator = sum(n), 
    .by = c(date, country)
  ) 

output <- data.frame()

for(c in unique(summarised$country)){
  
  just_country <- summarised |>
    dplyr::filter(country == c)
  
  min_date <- just_country |>
    dplyr::mutate(month = lubridate::floor_date(date, 'month')) |>
    dplyr::mutate(total = sum(numerator), .by = 'month') |>
    dplyr::filter(numerator >= 1 & total >= 2) |>
    dplyr::filter(date == min(date)) |>
    dplyr::pull(date)
  
  to_model <- just_country |>
    dplyr::filter(date >= min_date)
  
  output <- growth_rate(
    numerator = to_model$numerator,
    denominator = to_model$denominator,
    date = to_model$date,
    k_scaling = 10
    ) |>
    dplyr::left_join(just_country) |>
    dplyr::bind_rows(output)
  
}

#########################
# Visualize
########################

vec_for_plot <- c(-2, -5, -7, -10, -14, -21, -28, Inf, 28, 21, 14, 10, 7, 5, 2)

cols <- c('#007C91', '#582C83')

panel_A <- output |> 
  ggplot2::ggplot(
    ggplot2::aes(x = date, col = country, fill = country)
  )+
  ggplot2::geom_point(
    ggplot2::aes(y = numerator / denominator), pch = 21, col = 'black', show.legend = F, alpha = 0.5
  )+
  ggdist::geom_lineribbon(
    ggplot2::aes(y = prevalence, ymin = prevalence_lower, ymax = prevalence_upper), alpha = 0.75, size = 0.5, col = 'black'
  )+
  ggplot2::scale_fill_manual(
    name = '',
    values = cols
    )+
  ggplot2::scale_color_manual(
    name = '',
    values = cols
    )+
  ggplot2::scale_x_date(
    name = "",
    labels = scales::label_date_short(),
    date_breaks = "1 month"
  ) +
  ggplot2::scale_y_continuous(
    name = 'Variant share (proportion)'
    )+
  ggplot2::labs(title = glue::glue('{lineage} variant share
                                   
                                   '))+
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))

panel_B <- output |> 
  ggplot2::ggplot(
    ggplot2::aes(x = date, fill = country, col = country)
  )+
  ggplot2::geom_hline(
    ggplot2::aes(yintercept = 0), lty = 2, size = 0.5
  ) + 
  ggdist::geom_lineribbon(
    ggplot2::aes(y = gr, ymin = gr_lower, ymax = gr_upper), alpha = 0.75, size = 0.5
  )+
  ggplot2::scale_fill_manual(
    name = "",
    values = cols
  )+
  ggplot2::scale_color_manual(
    name = "",
    values = cols
  )+
  ggplot2::scale_x_date(
    name = "",
    labels = scales::label_date_short(),
    date_breaks = "1 month"
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))+
  ggplot2::scale_y_continuous(name = 'Daily relative growth rate', 
                              sec.axis = ggplot2::sec_axis(~., name="Doubling time", 
                                                           labels = function(x) log(2) / (x), 
                                                           breaks = log(2) / (vec_for_plot)))+
  ggplot2::labs(title = glue::glue('{lineage} growth rates
                                   
                                   '))
  
ggplot2::ggsave(panel_A + panel_B, file = glue::glue('Collab/outputs/comparison_figure_{lineage}.png'), dpi = 300, width = 12, height = 6)
