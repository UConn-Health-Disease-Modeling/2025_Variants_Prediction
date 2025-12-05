
###############################
# UKHSA / KCDC / UConn collab #
###############################
# Modelling relative growth rates

folder_wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)

# install packages
# source(glue::glue('{folder_wd}/variant_modelling/functions/setup.R'))

source(glue::glue('{folder_wd}/variant_modelling/functions/functions.R'))

# summary pull from GISAID: https://gisaid.org/hcov19-variants/
dataset <- readr::read_csv(glue::glue('{folder_wd}/variant_modelling/data/summary_gisaid.csv'))

# Each row is a country , day, lineage, and count of total tests sequenced (denominator) and those for that lineage (numerator)
head(dataset)

# Test run - compare lineages for SK and UK
############################################

# SETUP

# Pick a lineage: NB - code it setup to take the aliased variant name, as this is more common
modelled_lineage <- 'EG.5.1'

# the smoothness of the generalized additive model smoother
# smaller numbers are less smooth, greater numbers more smooth
k_scaling <- 14

# Optional arguments to set plot limits (useful for lineages that have been around for ages)
# defaults to no date limits for plotting
plot_lower_lim <- min(dataset$date)
plot_upper_lim <- max(dataset$date)

# axis increments
axis_increments <- '1 month'

# colours for each country
cols <- c('#007C91', '#582C83')

# Run the modelling

# As Pangolin assigns an alias to longer lineage names, I have included a function that
# gives you the alias for any lineage. This is useful for identifying sublineages that
# may end up with a second alias

first_part <- gsub(
  pattern = "[^a-zA-Z]",
  replacement = "",
  x = modelled_lineage
)

second_part <- gsub(
  pattern = "[a-zA-Z]",
  replacement = "",
  x = modelled_lineage
)

alias_list <- get_alias()

full_lineage <- alias_list |>
  dplyr::filter(
  alias == first_part
) |>
  dplyr::mutate(full_lineage = paste0(alias_lineage, second_part)) |>
  dplyr::pull(full_lineage)

if(length(full_lineage) == 0) full_lineage <- modelled_lineage

full_lineage

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
    is_target_lineage = stringr::str_starts(unaliased_lineage, paste0(full_lineage, '.')) | unaliased_lineage == full_lineage
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

  # model from the point where there are at least two sequenced cases in a week as
  # somtimes dates are incorrect
  min_date <- just_country |>
    dplyr::mutate(
      month = lubridate::floor_date(date, 'month')
      ) |>
    dplyr::mutate(
      total = sum(numerator),
      .by = 'month'
      ) |>
    dplyr::filter(
      numerator >= 1 & total >= 2
      ) |>
    dplyr::filter(
      date == min(date)
      ) |>
    dplyr::pull(
      date
      )

  to_model <- just_country |>
    dplyr::filter(
      date >= min_date
      )

  output <- growth_rate(
    numerator = to_model$numerator,
    denominator = to_model$denominator,
    date = to_model$date,
    k_scaling = k_scaling
  ) |>
    dplyr::left_join(just_country) |>
    dplyr::bind_rows(output)

}

#########################
# Visualize
########################

# prev and model fit, plot 1

panel_A <- output |>
  dplyr::filter(date >= plot_lower_lim, date <= plot_upper_lim) |>
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
    date_breaks = axis_increments
  ) +
  ggplot2::scale_y_continuous(
    name = 'Variant share (proportion)'
  )+
  ggplot2::labs(title = glue::glue('{modelled_lineage} variant share

                                   '))+
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))

# Growth rate plot
# as log(2) / GR gives us the doubling time, we can annotate this on the LHS axis

vec_for_plot <- c(-2, -5, -7, -10, -14, -21, -28, Inf, 28, 21, 14, 10, 7, 5, 2)

panel_B <- output |>
  dplyr::filter(date >= plot_lower_lim, date <= plot_upper_lim) |>
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
    date_breaks = axis_increments
  ) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))+
  ggplot2::scale_y_continuous(name = 'Daily relative growth rate',
                              sec.axis = ggplot2::sec_axis(~., name="Doubling time",
                                                           labels = function(x) log(2) / x,
                                                           breaks = log(2) / vec_for_plot))+
  ggplot2::labs(title = glue::glue('{modelled_lineage} growth rates

                                   '))

ggplot2::ggsave(panel_A + panel_B, file = glue::glue('{folder_wd}/variant_modelling/outputs/comparison_figure_{modelled_lineage}.png'), dpi = 300, width = 12, height = 6)
