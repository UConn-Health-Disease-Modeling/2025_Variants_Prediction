#!/usr/bin/env Rscript

# Draw stacked variant-share trajectories for all 15 selected countries.
#
# Run from the code directory:
#   cd code
#   Rscript draw_variants.R
#
# For each country, variants are ranked by the sum of their daily shares over
# the displayed period. Only the top 15 are coloured; "Others" is the summed
# share of every remaining variant in the input RDS. It is not calculated as
# one minus the selected variants. Each three-country figure uses one colour
# mapping, so the same variant has the same colour within that figure.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

figure_groups <- list(
  US_UK_SK = c("United States", "United Kingdom", "South Korea"),
  CA_DE_FR = c("Canada", "Germany", "France"),
  DK_AU_JP = c("Denmark", "Australia", "Japan"),
  ES_NL_IT = c("Spain", "Netherlands", "Italy"),
  SE_CH_BE = c("Sweden", "Switzerland", "Belgium")
)
countries <- unname(unlist(figure_groups, use.names = FALSE))
top_n <- 15L
smoothing_days <- 5L
plot_start <- as.IDate("2021-01-01")

# who2.rds is the historical name. variants_who2.rds and variants.rds are
# accepted so the script also works with the current simplified pipeline.
input_candidates <- c("who2.rds", "variants_who2.rds", "variants.rds")
available_inputs <- input_candidates[file.exists(input_candidates)]
if (length(available_inputs) == 0L) {
  stop(
    "Missing input. Expected one of: ",
    paste(input_candidates, collapse = ", ")
  )
}
input_rds <- available_inputs[1L]

mapping_file <- file.path(
  "..", "result", "02_group_variants", "global_lineage_mapping.csv"
)

rolling_mean_partial <- function(x, width) {
  left <- floor((width - 1L) / 2L)
  right <- width - left - 1L

  vapply(seq_along(x), function(i) {
    mean(x[max(1L, i - left):min(length(x), i + right)], na.rm = TRUE)
  }, numeric(1L))
}

# Draw ordinary legend keys without borders and the transparent excluded key
# with the same 0.25-mm black border used around the excluded plot region.
draw_variant_legend_key <- function(data, params, size) {
  key_fill <- if (!is.null(data$fill) && length(data$fill) > 0L) {
    data$fill[1L]
  } else {
    "transparent"
  }
  alpha_value <- grDevices::col2rgb(key_fill, alpha = TRUE)[4L, 1L]

  grid::rectGrob(
    gp = grid::gpar(
      col = if (alpha_value == 0L) "#222222" else NA_character_,
      fill = key_fill,
      lwd = 0.25 * 72.27 / 25.4
    )
  )
}

variant_data <- as.data.table(readRDS(input_rds))
required_columns <- c("country", "date", "variant", "numerator", "share")
missing_columns <- setdiff(required_columns, names(variant_data))
if (length(missing_columns) > 0L) {
  stop(
    input_rds, " is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

variant_data[, date := as.IDate(date)]
variant_data <- variant_data[
  country %chin% countries &
    date >= plot_start &
    !is.na(variant) &
    tolower(trimws(variant)) != "unassigned"
]

missing_countries <- setdiff(countries, unique(variant_data$country))
if (length(missing_countries) > 0L) {
  stop("Missing countries in ", input_rds, ": ", paste(missing_countries, collapse = ", "))
}

plot_end <- max(variant_data$date, na.rm = TRUE)

daily_total_check <- variant_data[, .(
  total_share = sum(share, na.rm = TRUE)
), by = .(country, date)][total_share > 1 + 1e-10]
if (nrow(daily_total_check) > 0L) {
  stop("The summed input-RDS share exceeds 100% for at least one country-date.")
}

# Select each country's top 15 by cumulative daily share. Sequence counts and
# stable tie-breakers make the ranking deterministic when sums are equal.
variant_totals <- variant_data[, .(
  sum_share = sum(share, na.rm = TRUE),
  sum_sequences = sum(numerator, na.rm = TRUE)
), by = .(country, variant)]
setorder(variant_totals, country, -sum_share, -sum_sequences, variant)
variant_totals[, rank := seq_len(.N), by = country]
top15 <- variant_totals[rank <= top_n]

top15_check <- top15[, .N, by = country]
if (any(top15_check$N != top_n)) {
  stop("Every country must have at least ", top_n, " variants.")
}

# Use concise Pangolin aliases in the legend when the current pipeline mapping
# is available. The underlying variant value remains the colour key.
display_mapping <- unique(top15[, .(
  variant,
  display_label = variant
)])

mapping_applied <- FALSE

if (file.exists(mapping_file)) {
  lineage_mapping <- fread(mapping_file)
  mapping_columns <- c("original_lineage", "canonical_lineage", "harmonized_group")

  if (all(mapping_columns %in% names(lineage_mapping))) {
    lineage_mapping[, exact_group_label := canonical_lineage == harmonized_group]
    lineage_mapping[, label_length := nchar(original_lineage)]
    setorder(
      lineage_mapping,
      harmonized_group,
      -exact_group_label,
      label_length,
      original_lineage
    )

    short_labels <- lineage_mapping[
      nzchar(original_lineage),
      .SD[1L],
      by = harmonized_group
    ][, .(
      variant = harmonized_group,
      short_label = original_lineage
    )]

    display_mapping[short_labels, on = "variant", display_label := i.short_label]
    mapping_applied <- TRUE
  }
}

# If the diagnostic mapping CSV is unavailable, recover concise aliases from
# the same frozen Pango reference used by 02_select_variants.R. This keeps the
# figure labels stable without changing the harmonized values used for colour
# assignment or aggregation.
if (!mapping_applied) {
  pangolin_bin <- Sys.getenv(
    "PANGOLIN_BIN",
    unname(Sys.which("pangolin"))
  )

  if (nzchar(pangolin_bin) &&
      file.exists(pangolin_bin) &&
      requireNamespace("jsonlite", quietly = TRUE)) {
    alias_text <- system2(pangolin_bin, "--aliases", stdout = TRUE, stderr = FALSE)
    alias_definitions <- jsonlite::fromJSON(
      paste(alias_text, collapse = "\n"),
      simplifyVector = FALSE
    )

    expand_alias <- function(lineage, max_steps = 30L) {
      current <- lineage
      for (step in seq_len(max_steps)) {
        prefix <- sub("\\..*$", "", current)
        if (prefix %chin% c("A", "B")) return(current)
        definition <- alias_definitions[[prefix]]
        if (is.null(definition) ||
            !is.character(definition) ||
            length(definition) != 1L ||
            !nzchar(definition)) {
          return(current)
        }
        current <- paste0(
          definition,
          substring(current, nchar(prefix) + 1L)
        )
      }
      current
    }

    reversible_aliases <- setdiff(names(alias_definitions), c("A", "B"))
    reversible_aliases <- reversible_aliases[vapply(
      alias_definitions[reversible_aliases],
      function(x) is.character(x) && length(x) == 1L && nzchar(x),
      logical(1L)
    )]
    expanded_alias_roots <- vapply(
      reversible_aliases,
      expand_alias,
      character(1L)
    )

    concise_alias <- function(lineage) {
      matches <- which(
        lineage == expanded_alias_roots |
          startsWith(lineage, paste0(expanded_alias_roots, "."))
      )
      if (length(matches) == 0L) return(lineage)

      candidates <- paste0(
        reversible_aliases[matches],
        substring(lineage, nchar(expanded_alias_roots[matches]) + 1L)
      )
      candidates <- c(lineage, candidates)
      candidates[order(nchar(candidates), candidates)][1L]
    }

    raw_lineage_file <- file.path(
      "..",
      "reference",
      "UKHSA-UConn-variant-modelling",
      "variant_modelling",
      "data",
      "summary_GISAID_20240918.csv"
    )

    if (file.exists(raw_lineage_file)) {
      observed_labels <- unique(fread(
        raw_lineage_file,
        select = "lineage"
      ))
      observed_labels <- observed_labels[
        !is.na(lineage) &
          nzchar(trimws(lineage)) &
          tolower(trimws(lineage)) != "unassigned"
      ]
      observed_labels[, canonical_lineage := vapply(
        lineage,
        expand_alias,
        character(1L)
      )]

      observed_display_label <- function(group) {
        candidates <- observed_labels[
          canonical_lineage == group |
            startsWith(canonical_lineage, paste0(group, "."))
        ]
        if (nrow(candidates) == 0L) return(concise_alias(group))

        candidates[, exact_group := canonical_lineage == group]
        candidates[, label_length := nchar(lineage)]
        setorder(candidates, -exact_group, label_length, lineage)
        candidates$lineage[1L]
      }

      display_mapping[, display_label := vapply(
        variant,
        observed_display_label,
        character(1L)
      )]
    } else {
      display_mapping[, display_label := vapply(
        variant,
        concise_alias,
        character(1L)
      )]
    }
  }
}

# Prevent two different variants from receiving indistinguishable legend text.
display_mapping[, display_label := ifelse(
  duplicated(display_label) | duplicated(display_label, fromLast = TRUE),
  paste0(display_label, " (", variant, ")"),
  display_label
)]

top15[display_mapping, on = "variant", display_label := i.display_label]
top15[, country_order := match(country, countries)]
setorder(top15, country_order, rank)
top15[, country_order := NULL]

# Keep the non-top-15 variants as the original grey "Others" band. Only the
# unrepresented space above that band, through 100%, is transparent.
others_colour <- "#D9D9D9"

make_figure_palette <- function(figure_countries) {
  selected_variants <- sort(unique(
    top15[country %chin% figure_countries]$variant
  ))

  setNames(
    grDevices::hcl(
      h = ((seq_along(selected_variants) - 1L) * 137.508) %% 360,
      c = 82,
      l = rep(c(50, 66), length.out = length(selected_variants)),
      fixup = TRUE
    ),
    selected_variants
  )
}

make_country_plot_data <- function(country_name) {
  selected <- copy(top15[country == country_name][order(rank)])
  country_variants <- selected$variant
  calendar <- data.table(date = seq(plot_start, plot_end, by = "day"))

  daily <- variant_data[
    country == country_name & variant %chin% country_variants,
    .(share = sum(share, na.rm = TRUE)),
    by = .(date, variant)
  ]

  plot_data <- CJ(
    date = calendar$date,
    variant = country_variants,
    unique = TRUE
  )
  plot_data[daily, on = .(date, variant), share := i.share]
  plot_data[is.na(share), share := 0]

  # "Others" contains only the remaining variants present in who2/current
  # variants.rds. Share absent from that input stays blank rather than being
  # filled up to 100%.
  others_daily <- variant_data[
    country == country_name & !variant %chin% country_variants,
    .(share = sum(share, na.rm = TRUE)),
    by = date
  ]
  others <- copy(calendar)
  others[others_daily, on = "date", share := i.share]
  others[is.na(share), share := 0]
  others[, variant := "Others"]

  plot_data <- rbindlist(list(plot_data, others), use.names = TRUE)
  setorder(plot_data, variant, date)
  plot_data[, share_smooth := rolling_mean_partial(share, smoothing_days), by = variant]
  plot_data[!is.finite(share_smooth), share_smooth := 0]

  # Preserve the absolute share represented in the input RDS. Do not normalize
  # the coloured and grey bands to 100%.
  plot_data[, share_plot := pmax(0, share_smooth)]

  peak_order <- plot_data[
    variant != "Others",
    .(peak_date = date[which.max(share_plot)][1L]),
    by = variant
  ][order(peak_date, variant)]$variant
  stack_order <- c(peak_order, "Others")

  plot_data[, stack_rank := match(variant, stack_order)]
  setorder(plot_data, date, stack_rank)
  plot_data[, ymax := pmin(1, pmax(0, cumsum(share_plot))), by = date]
  plot_data[, ymin := pmin(1, pmax(0, ymax - share_plot))]

  list(data = plot_data, selected = selected)
}

make_country_plot <- function(country_name, variant_colours) {
  prepared <- make_country_plot_data(country_name)
  excluded_label <- "Excluded(D)"
  legend_variants <- c(prepared$selected$variant, "Others", excluded_label)
  legend_labels <- c(prepared$selected$display_label, "Others", excluded_label)
  country_colours <- c(
    variant_colours[prepared$selected$variant],
    Others = others_colour,
    setNames("transparent", excluded_label)
  )

  # The transparent region begins at the upper boundary of the grey "Others"
  # band and extends to 100%, with a thin black outline.
  transparent_region <- prepared$data[
    variant == "Others",
    .(date, ymin = ymax, ymax = 1)
  ]
  transparent_verticals <- transparent_region[
    c(1L, nrow(transparent_region))
  ]
  transparent_top <- data.table(
    date = c(
      transparent_region$date[1L],
      transparent_region$date[nrow(transparent_region)]
    ),
    y = 1
  )

  ggplot(
    prepared$data,
    aes(
      x = as.Date(date),
      ymin = ymin,
      ymax = ymax,
      fill = variant,
      group = variant
    )
  ) +
    geom_ribbon(colour = NA, key_glyph = draw_variant_legend_key) +
    geom_line(
      data = transparent_region,
      aes(x = as.Date(date), y = ymin),
      inherit.aes = FALSE,
      colour = "#222222",
      linewidth = 0.25,
      lineend = "butt"
    ) +
    geom_line(
      data = transparent_top,
      aes(x = as.Date(date), y = y),
      inherit.aes = FALSE,
      colour = "#222222",
      linewidth = 0.25,
      lineend = "butt"
    ) +
    geom_segment(
      data = transparent_verticals,
      aes(
        x = as.Date(date),
        xend = as.Date(date),
        y = ymin,
        yend = ymax
      ),
      inherit.aes = FALSE,
      colour = "#222222",
      linewidth = 0.25,
      lineend = "butt"
    ) +
    scale_fill_manual(
      values = country_colours,
      breaks = legend_variants,
      labels = legend_labels,
      limits = legend_variants,
      drop = FALSE,
      guide = guide_legend(title = NULL, ncol = 1, byrow = TRUE)
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_x_date(
      limits = as.Date(c(plot_start, plot_end)),
      date_breaks = "1 year",
      date_labels = "%Y",
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(clip = "off") +
    labs(title = country_name, x = NULL, y = "Domestic share") +
    theme_classic(base_size = 13) +
    theme(
      axis.line = element_line(linewidth = 0.45, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#222222"),
      axis.text = element_text(colour = "#222222", size = 10),
      axis.title.y = element_text(size = 11, margin = margin(r = 8)),
      plot.title = element_text(size = 14, face = "bold", hjust = 0),
      legend.text = element_text(size = 9),
      legend.key.height = grid::unit(4.0, "mm"),
      legend.key.width = grid::unit(4.0, "mm"),
      legend.box.margin = margin(0, 0, 0, 8),
      plot.margin = margin(7, 12, 7, 7)
    )
}

make_figure <- function(group_name, figure_countries) {
  variant_colours <- make_figure_palette(figure_countries)
  combined_plot <- wrap_plots(
    lapply(
      figure_countries,
      make_country_plot,
      variant_colours = variant_colours
    ),
    ncol = 1,
    heights = rep(1, length(figure_countries))
  )
  output_pdf <- paste0(group_name, ".pdf")

  ggsave(
    output_pdf,
    combined_plot,
    width = 16,
    height = 10.5,
    units = "in",
    device = grDevices::pdf,
    useDingbats = FALSE,
    bg = "transparent"
  )

  output_pdf
}

output_files <- vapply(
  names(figure_groups),
  function(group_name) make_figure(group_name, figure_groups[[group_name]]),
  character(1L)
)

message("Input: ", input_rds)
message("Saved figures: ", paste(unname(output_files), collapse = ", "))
