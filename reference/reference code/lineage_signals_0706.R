
###############################
# UKHSA / KCDC / UConn collab #
###############################
# Find groups of lineages to model

# folder_wd <- "C://Users//young//OneDrive//Desktop//UKHSA_COVID19_Variants"
#folder_wd <- system("echo $(git rev-parse --show-toplevel)/", intern = TRUE)


# install packages
# source(glue::glue('{folder_wd}/variant_modelling/functions/setup.R'))

source(glue::glue('/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants/Code/functions_0706.R'))

alias_list <- get_alias()

# summary pull from GISAID: https://gisaid.org/hcov19-variants/
# Pull out just one country for this assessment
# Remove unassigned lineages
# Set a date limit on the total time series
country_select <- 'United Kingdom'

months_from_most_recent <- 6* 4 * 7

dataset <- readr::read_csv(glue::glue('/Users/frankyzhang/Dropbox/Jo_Franky/2024_Variants/Code/summary_gisaid_new.csv')) |>
  dplyr::filter(country == country_select,
                lineage != 'Unassigned',
                date >= max(date, na.rm = T) - months_from_most_recent) |>
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


# What proportion should a lineage be represented at before we model it
p_lim <- 0.015
# How many sequenced cases do we require before we model a lineage group
n_lim <- 50
# over what range of time from the most recent date should this be applied?
date_range <- 6 * 7

summarised_two_weeks <- dataset |>
  dplyr::filter(date >= (max(date) - date_range)) |>
  dplyr::summarise(n = sum(n), .by = unaliased_lineage)

# How many of these signals in the previous weeks should we monitor?

classification_to_append <- get_classification(
  lineages = summarised_two_weeks$unaliased_lineage,
  number_sequences = summarised_two_weeks$n,
  p_lim = p_lim, n_lim = n_lim, alias_list
) |>
  dplyr::mutate(
    decimal_lineage = paste0(classified_unasliased, '.')
  )

classified_data <- dataset |>
  dplyr::mutate(
    decimal_lineage = paste0(unaliased_lineage, '.')
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
    length_lineage = stringr::str_count(unaliased_lineage, '\\.'),
    classified_label = dplyr::if_else((length_lineage - length_class) <= 1,
                                      classified_label,
                                      NA),
    classified_label = dplyr::if_else(is.na(classified_label),
                                      "Other",
                                      classified_label)
  ) |>
  dplyr::filter(classified_label != 'Other')

lineages_to_model <- unique(classified_data$classified_label)

output <- data.frame()

for(lin in lineages_to_model){

  this_lineage <- classified_data |>
    dplyr::mutate(denominator = sum(n), .by = date) |>
    dplyr::summarise(numerator = sum(n * (classified_label == lin)),
                     .by = c(date, denominator))

  min_date <- this_lineage |>
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

  to_model <- this_lineage |>
    dplyr::filter(
      date >= min_date
    )

  output <- growth_rate(
    numerator = to_model$numerator,
    denominator = to_model$denominator,
    date = to_model$date,
    k_scaling = 10
  ) |>
    dplyr::mutate(classified_label = lin) |>
    dplyr::left_join(this_lineage) |>
    dplyr::bind_rows(output)


}

p1 <- output |>
  ggplot2::ggplot(
    ggplot2::aes(x = date)
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = numerator / denominator), pch = 21, col = 'black', show.legend = F, alpha = 0.5, fill = '#007C91'
  )+
  ggdist::geom_lineribbon(
    ggplot2::aes(y = prevalence, ymin = prevalence_lower, ymax = prevalence_upper), alpha = 0.75, size = 0.5, col = 'black', fill = '#007C91'
  )+
  ggplot2::facet_wrap(
    ~ classified_label
  )+
  ggplot2::scale_x_date(
    name = "",
    labels = scales::label_date_short(),
    date_breaks = '4 weeks'
  ) +
  ggplot2::scale_y_continuous(
    name = 'Lineage share (proportion)'
  )+
  ggplot2::labs(title = glue::glue('lineage share'))

vec_for_plot <- c(-2, -5, -7, -10, -14, -21, Inf, 21, 14, 10, 7, 5, 2)


p2 <- output |>
  ggplot2::ggplot(
    ggplot2::aes(x = date)
  )+
  ggplot2::geom_hline(
    ggplot2::aes(yintercept = 0), lty = 2, size = 0.5
  ) +
  ggdist::geom_lineribbon(
    ggplot2::aes(y = gr, ymin = gr_lower, ymax = gr_upper), alpha = 0.75, size = 0.5, fill = '#007C91'
  )+
  ggplot2::scale_x_date(
    name = "",
    labels = scales::label_date_short(),
    date_breaks = '4 weeks'
  ) +
  ggplot2::facet_wrap(
    ~ classified_label
  )+
  ggplot2::scale_y_continuous(name = 'Daily relative growth rate',
                              sec.axis = ggplot2::sec_axis(~., name="Doubling time",
                                                           labels = function(x) log(2) / x,
                                                           breaks = log(2) / vec_for_plot))+
  ggplot2::labs(title = glue::glue('lineage growth rates

                                   '))

final_date <- output |>
  dplyr::filter(!is.na(gr)) |>
  dplyr::filter(date == max(date))

p3 <- final_date |>
  ggplot2::ggplot(
    ggplot2::aes(y = forcats::fct_reorder(classified_label, -prevalence))
  )+

  ggplot2::geom_errorbarh(
    ggplot2::aes(x = prevalence, xmin = prevalence_lower , xmax = prevalence_upper           ),  height = 0.2
  )+
  ggplot2::geom_point(ggplot2::aes(x = prevalence), pch = 21, size = 3, col = 'black', show.legend = F, fill = '#007C91')+
  ggplot2::scale_y_discrete(
    name = ""
  ) +
  ggplot2::scale_x_continuous(
    name = 'Lineage share (proportion)'
  )+
  ggplot2::labs(title = glue::glue('{country_select} lineage share on {unique(final_date$date)}

                                   '))

p4 <- final_date |>
  ggplot2::ggplot(
    ggplot2::aes(y = forcats::fct_reorder(classified_label, -gr))
  )+
  ggplot2::geom_vline(
    ggplot2::aes(xintercept = 0), lty = 2, size = 0.5
  ) +
  ggplot2::geom_errorbarh(
    ggplot2::aes(x = gr, xmin = gr_lower, xmax = gr_upper),  height = 0.2
  )+
  ggplot2::geom_point(ggplot2::aes(x = gr), pch = 21, size = 3, col = 'black', show.legend = F, fill = '#007C91')+
  ggplot2::scale_y_discrete(
    name = ""
  ) +
  ggplot2::scale_x_continuous(name = 'Daily relative growth rate')+
  ggplot2::labs(title = glue::glue('{country_select} growth rates on {unique(final_date$date)}

                                   '))


combined_plot <- p1 + p2 + p3 + p4 + plot_layout(design = c('AAAC
                                                             BBBD'))

ggplot2::ggsave(combined_plot, file = glue::glue('C://Users//young//OneDrive//Desktop//UKHSA_COVID19_Variants//outputs//{country_select}_lineages_{unique(final_date$date)}.png'),
                dpi = 300, width = 18, height = 16)

write.csv(classified_data, "C://Users//young//OneDrive//Desktop//UKHSA_COVID19_Variants//classified_data_uk.csv")
