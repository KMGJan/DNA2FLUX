#!/usr/bin/env Rscript

# Thu Jun 26 10:28:10 2025 ------------------------------

# This scripts contains all the helper functions

# Import files
readAndArrange <- function(file = file, arrange = TRUE) {
  df <- read_csv(file = file, show_col_types = FALSE)
  if (arrange) df |> arrange(sample_week) else df
}

# Assign season based on ISO week
add_season <- function(df) {
  # Spring: week 5-25, Summer: 26-37, Fall: 38-50, Winter the remaining weeks
  df |>
    mutate(
      iso_week = isoweek(sample_week),
      season = case_when(
        iso_week %in% 1:10 ~ "Winter",
        iso_week %in% 11:22 ~ "Spring",
        iso_week %in% 23:35 ~ "Summer",
        iso_week %in% 36:48 ~ "Fall",
        iso_week %in% 49:53 ~ "Winter",
        TRUE ~ "Winter"
      )
    )
}

# Ensure season is a factor with correct order
factor_season <- function(df) {
  df |>
    mutate(
      season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))
    )
}

# Filter by iso week and year, optionally summarise biomass
filter_my_df <- function(df, summarise = TRUE) {
  df <- df |>
    filter(isoweek(sample_week) %in% 2:51, year(sample_week) > 2007)
  if (summarise) {
    df |>
      group_by(node_name, season, station_name) |>
      summarise(avg = mean(biomass, na.rm = TRUE), .groups = "drop")
  } else {
    df
  }
}

# PCoA helper
process_pcoa <- function(pcoa_result, metadata, id_cols = 1:3) {
  coords <- as.data.frame(pcoa_result$vectors[, 1:2])
  names(coords) <- c("Axis1", "Axis2")

  coords |>
    bind_cols(metadata[, id_cols]) |>
    mutate(
      iso_week = isoweek(sample_week),
      season = case_when(
        iso_week %in% 1:10 ~ "Winter",
        iso_week %in% 11:22 ~ "Spring",
        iso_week %in% 23:35 ~ "Summer",
        iso_week %in% 36:48 ~ "Autumn",
        iso_week %in% 49:53 ~ "Winter",
        TRUE ~ NA_character_
      ),
      month = month(sample_week),
      month_abb = factor(month.abb[month], levels = month.abb)
    )
}
PCOA_plot <- function(df, eig) {
  df |>
    ggplot(aes(x = Axis1, y = Axis2)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_point(size = 2, mapping = aes(fill = predator, shape = predator)) +
    facet_wrap(~month_abb) +
    scale_shape_manual(values = c(21, 22, 24)) +
    scale_fill_manual(values = color_mapping) +
    coord_equal(xlim = c(-0.7, 0.7), ylim = c(-0.6, 0.6)) +
    scale_x_continuous(breaks = seq(-1, 1, 0.5)) +
    scale_y_continuous(breaks = seq(-1, 1, 0.5)) +
    labs(
      x = paste0("Axis 1 (", round(100 * eig[1], 1), "%)"),
      y = paste0("Axis 2 (", round(100 * eig[2], 1), "%)")
    )
}

# Helper to bind PCOA data with type
bind_pcoa <- function(pp, zp) {
  bind_rows(
    pp |> mutate(type = "phytoplankton"),
    zp |> mutate(type = "zooplankton")
  )
}


# Diet overlap
interspecific_overlap <- function(mat, cols, method = "schoener") {
  # Create numeric matrix and row names
  mat2 <- mat[, cols] |> as.matrix()
  rownames(mat2) <- with(mat, paste(station, predator, sample_week, sep = "_"))

  # Get Schoener's D overlap matrix and long format
  niche_df <- niche.overlap(t(mat2), method = method) |>
    dist2list() |>
    as.data.frame() |>
    #mutate(value = ifelse(col == row, NA, as.numeric(value))) |>
    separate(
      col,
      into = c("station.x", "predator.x", "sample_week.x"),
      sep = "_"
    ) |>
    separate(
      row,
      into = c("station.y", "predator.y", "sample_week.y"),
      sep = "_"
    ) |>
    mutate(
      x = predator.x,
      y = predator.y,
      xy_min = pmin(x, y),
      xy_max = pmax(x, y)
    ) |>
    filter(x == xy_min, y == xy_max) |>
    select(-xy_min, -xy_max) |>
    na.omit()
}

plot_interspecific_overlap <- function(data, ...) {
  data |>
    mutate(interaction = paste(x, y, sep = "_")) |>
    group_by(iso_week = isoweek(sample_week.x), interaction, facet, type) |>
    summarise(
      y = mean(value),
      ymin = pmax(0, y - sd(value)),
      ymax = pmin(1, y + sd(value)),
      .groups = "drop"
    ) |>
    ggplot(aes(
      x = iso_week,
      y = y,
      ymin = ymin,
      ymax = ymax,
      col = interaction,
      fill = interaction,
      group = interaction
    )) +
    geom_ribbon(alpha = .2) +
    geom_line() +

    labs(x = NULL, y = "Schoener's D") +
    scale_x_continuous(
      breaks = seq(1, 52.1775, 4.348125),
      labels = month.abb,
      expand = c(0, 0)
    ) +
    facet_grid(. ~ facet)
}


# Reusable plot function with and without forage ratio
plot_flux_difference <- function(trophic_level_filter) {
  isoweek_fluxes |>
    left_join(
      isoweek_null |> rename("null" = avg_fluxes),
      by = c("iso_week", "station", "predator", "prey")
    ) |>
    mutate(
      predator = factor(predator, levels = node_data$node_name),
      prey = factor(prey, levels = node_data$node_name),
      trophic_level = case_when(
        prey %in% node_data$node_name[node_data$trophic_level == 2] ~
          "zooplankton",
        prey %in% node_data$node_name[node_data$trophic_level == 1] ~
          "phytoplankton"
      )
    ) |>
    filter(trophic_level == trophic_level_filter, iso_week %in% 2:51) |>
    na.omit() |>
    ggplot(aes(x = iso_week)) +
    geom_ribbon(
      aes(y = mean, ymin = lower, ymax = upper),
      alpha = .2,
      linewidth = .2,
      fill = "black",
      col = "black"
    ) +
    geom_line(aes(y = mean, color = "Selectivity")) +
    geom_line(aes(y = null, color = "No selectivity")) +
    geom_line(aes(y = mean, color = "Selectivity")) +
    facet_grid(
      prey ~ predator,
      scales = "free_y"
    ) +
    scale_color_manual(
      name = NULL,
      guide = "none",
      values = c("Selectivity" = "black", "No selectivity" = "#ff7f00")
    ) +
    scale_x_continuous(
      breaks = seq(1, 52.1775, 4.348125 * 2),
      labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"),
      expand = c(0, 0)
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    labs(
      y = "Fluxes \n [kJ/day/m2]",
      x = NULL
    )
}
# Over the entire timeseries
# # ==== Plotting function ====
plot_flux_diff <- function(data, level, ylim) {
  df <- data |> filter(trophic_level == level)

  prey_levels <- levels(droplevels(df$prey))
  bg_df <- tibble(
    prey = factor(prey_levels, levels = prey_levels),
    idx = seq_along(prey_levels)
  ) |>
    filter(idx %% 2 == 0)

  dodge <- position_dodge2(width = 0.4, preserve = "single")

  ggplot(
    df,
    aes(
      x = prey,
      y = rel_mean,
      ymax = rel_upper,
      ymin = rel_lower,
      fill = predator,
      alpha = sig
    )
  ) +
    geom_rect(
      data = bg_df,
      aes(
        xmin = as.numeric(prey) - 0.5,
        xmax = as.numeric(prey) + 0.5,
        ymin = -Inf,
        ymax = Inf
      ),
      fill = "gray90",
      inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    geom_bar(stat = "identity", position = dodge) +
    geom_errorbar(position = dodge) +
    scale_fill_manual(values = color_mapping) +
    scale_alpha_manual(values = c(0.4, 1)) +
    scale_y_continuous(limits = ylim) +
    coord_flip() +
    labs(
      y = "Difference from null model [kJ/day/m2]",
      x = NULL
    )
}
# linear models
# Run the linear model for each trophic level using broom and purrr
run_linear_model_run <- function(df, group, formula_expr, ...) {
  formula_parsed <- as.formula(formula_expr)
  df |>
    group_by({{ group }}) |>
    nest() |>
    mutate(
      model = map(data, ~ lm(formula_parsed, data = .x)),
      tidied = map(model, tidy),
      glanced = map(model, glance)
    )
}

# Save the residual plots
save_residuals <- function(mod, group_col, path.end, plot.title, ...) {
  mod |>
    select(group = all_of(group_col), model) |>
    pwalk(function(group, model) {
      pdf(
        file = file.path(
          "output",
          "residuals",
          paste0(group, path.end)
        )
      )
      par(mfrow = c(2, 2))
      plot(model)
      mtext(
        paste(group, plot.title),
        outer = TRUE,
        line = -1.5,
        cex = 1.2
      )
      dev.off()
    })
}

# Summarise the model outputs
summarise_model <- function(mod, slope = "z", ...) {
  mod |>
    mutate(
      glanced = map(
        glanced,
        ~ select(.x, r.squared, p.value)
      ),
      tidied = map(
        tidied,
        ~ select(.x, term, estimate) |>
          mutate(term = if_else(term == slope, "slope", "intercept"))
      )
    ) |>
    unnest(c(glanced)) |>
    mutate(
      r.squared = round(r.squared, 3),
      P = if_else(p.value <= 0.05, "sig", "nosig")
    ) |>
    unnest(tidied)
}
# All workflow
run_model_workflow <- function(
  df,
  group_var,
  formula_expr,
  filename_suffix,
  plot_title,
  slope = "z"
) {
  mod <- run_linear_model_run(df, {{ group_var }}, formula_expr)

  save_residuals(
    mod = mod,
    group_col = as_label(enquo(group_var)),
    path.end = filename_suffix,
    plot.title = plot_title
  )

  summarise_model(mod, slope)
}
