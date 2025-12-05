
require(patchwork)

###############
# Plot theme #
##############

ggplot2::theme_set(
  ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                 axis.line.x = ggplot2::element_line(),
                 axis.line = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(face = 'bold'),
                 panel.grid.major = ggplot2::element_line(colour = 'white', linetype = 1, linewidth = 0.25),
                 panel.grid.minor = ggplot2::element_line(colour = 'white', linetype = 1, linewidth = 0.1),
                 panel.background =  ggplot2::element_rect(colour = 'grey95', fill ='grey95'),
                 legend.position = c(0.5, 0.99),
                 legend.background = ggplot2::element_blank()
  )

)

##################################
# Retrieve the latest alias list #
##################################

get_alias <- function() {

  myData <- rjson::fromJSON(file = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/pango_designation/alias_key.json")

  reference_list <- as.data.frame((myData[names(myData)[lapply(myData, length) <= 1]]))

  alias_list <- reference_list |>
    tidyr::gather(
      alias,
      alias_lineage
    ) |>
    dplyr::filter(
      !(alias_lineage %in% c("", "BA"))
    ) |>
    dplyr::mutate(
      alias_lineage = gsub(pattern = "B.1.1.529",
                           replacement = "BA",
                           alias_lineage)
    )

  return(alias_list)

}

###########################
# Inverse logit transform #
###########################

ilogit <- function(x) 1 / (1 + exp(- x))

########################
# Growth rate function #
########################

growth_rate <- function(
    numerator,
    denominator,
    date,
    k_scaling = NULL,
    basis = 'gp'
) {

  data <- tidyr::tibble(
    numerator = numerator,
    denominator = denominator,
    date = date
  ) |>
    dplyr::mutate(
      time_numeric = as.numeric(date - min(date)),
      time_factor = as.factor(time_numeric)
    )

  if (is.null(k_scaling)) {

    k <- ceiling(max(data$time_numeric) / 14)

  } else {

    k <- ceiling(max(data$time_numeric) / k_scaling)

  }

  if(basis == 'tp') "Thin plate splines ignore K, see `?s`"

  formula <- as.formula(
    cbind(numerator, denominator - numerator) ~
      s(time_numeric, k = k, bs = basis)
    + s(time_factor, bs = "re")
  )

  model <- mgcv::bam(
    formula = formula,
    data = data,
    discrete = T,
    family = "binomial"(
      link = "logit"
    )
  )

  ##################################
  # Model fit
  ##################################

  new_data <- expand.grid(
    date = seq(
      min(date),
      max(date),
      by = 1)
  ) |>
    dplyr::mutate(
      date = lubridate::ymd(date),
      time_numeric = as.numeric(date - min(date)),
      time_factor = as.factor(time_numeric)
    )

  fits <- mgcv::predict.gam(
    model,
    newdata = new_data,
    exclude = "s(time_factor)",
    se.fit = T
  )

  derivs <- gratia::derivatives(
    model,
    term = 's(time_numeric)',
    type = 'forward',
    eps = 1,
    n = nrow(new_data),
    new_data = new_data,
    frequentist = T
  ) |>
    dplyr::select(
      time_numeric = data,
      gr = derivative,
      gr_lower = lower,
      gr_upper = upper
    ) |>
    dplyr::slice(
      1,
      .by = time_numeric
    ) |>
    dplyr::mutate(
      time_numeric = floor(time_numeric)
    )

  new_data <- new_data |>
    dplyr::mutate(
      fit_logit = fits$fit,
      fit_logit_lower = fits$fit - 1.96 * fits$se.fit,
      fit_logit_upper = fits$fit + 1.96 * fits$se.fit,
      prevalence = ilogit(fit_logit),
      prevalence_lower = ilogit(fit_logit_lower),
      prevalence_upper = ilogit(fit_logit_upper)
    ) |>
    dplyr::left_join(
      derivs
    )

  to_return <- new_data |>
    dplyr::select(
      date,
      fit_logit,
      fit_logit_lower,
      fit_logit_upper,
      prevalence,
      prevalence_lower,
      prevalence_upper,
      gr,
      gr_lower,
      gr_upper
    )

  return(to_return)

}

##################################
# Group lineages based on
##################################

get_classification <- function(lineages, number_sequences, p_lim = 0.01, n_lim = 50, alias_list) {

  classifications <- tidyr::tibble(
    lineages = lineages,
    n = number_sequences,
    N = sum(number_sequences)
  ) |>
    dplyr::mutate(
      p = n / N
    ) |>
    dplyr::mutate(
      classified_interim = lineages,
      classified_final = dplyr::if_else(p > p_lim & n > n_lim,
                                        lineages,
                                        NA)
    )

  dim_old <- nrow(classifications)

  while (any(is.na(classifications$classified_final))) {

    classifications |>
      dplyr::pull(classified_interim) -> possible_parents

    found_in <- lapply(possible_parents, stringr::str_starts, pattern = paste0(possible_parents, "."))

    classifications |>
      dplyr::rowwise() |>
      dplyr::mutate(
        best_matches = list(possible_parents[stringr::str_starts(classified_interim,
                                                                 stringr::coll(paste0(possible_parents, ".")))])
      ) |>
      dplyr::mutate(
        selected_parent = ifelse(length(best_matches) == 0,
                                 NA,
                                 best_matches[which.max(stringr::str_length(best_matches))]),
        cropped = stringr::str_split(classified_interim,  stringr::coll(".")),
        cropped_processed = ifelse(length(cropped) > 2,
                                   paste((cropped)[-length((cropped))], collapse = "."),
                                   paste((cropped), collapse = ".")),
        new_list = ifelse(is.na(selected_parent),
                          cropped_processed,
                          selected_parent)
      )  |>
      dplyr::ungroup() |>
      dplyr::mutate(
        classified_interim = dplyr::if_else(is.na(classified_final),
                                            new_list,
                                            classified_final)
      ) |>
      dplyr::summarise(
        n = sum(n),
        .by = c("classified_interim", "N")
      ) |>
      dplyr::mutate(
        p = n / N
      ) |>
      dplyr::mutate(
        classified_final = dplyr::if_else(p > p_lim & n > n_lim,
                                          classified_interim,
                                          NA)
      ) -> classifications

    if (nrow(classifications) == dim_old) {

      break("Hunting over")

    }

    dim_old <- nrow(classifications)

  }

  classifications |>
    tidyr::drop_na() |>
    fuzzyjoin::fuzzy_left_join(
      alias_list,
      by = c("classified_final" = "alias_lineage"),
      match_fun = stringr::str_starts
    ) |>
    dplyr::mutate(
      alias_lineage = dplyr::if_else(is.na(alias_lineage),
                                     classified_final,
                                     alias_lineage),
      length_alias = stringr::str_length(alias_lineage)
    ) |>
    dplyr::filter(
      length_alias == max(length_alias),
      .by = classified_final
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      lineage_class_alias = gsub(x = classified_final,
                                 pattern = alias_lineage,
                                 replacement = "",
                                 fixed = T),
      classified_alias = dplyr::if_else(is.na(alias),#stringr::str_length(lineage_class_alias) == 0,
                                        classified_final,
                                        paste0(alias, lineage_class_alias)),
      classified_label = dplyr::if_else(is.na(alias),#stringr::str_length(lineage_class_alias) == 0,
                                        classified_final,
                                        paste0(classified_alias, " (", classified_final, ")"))
    ) |>
    dplyr::rename(
      classified_unasliased = classified_final
    ) |>
    dplyr::distinct(
      classified_unasliased,
      classified_alias,
      classified_label
    ) |>
    dplyr::ungroup() -> reference_frame

  return(reference_frame)

}

