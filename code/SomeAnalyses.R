#!/usr/bin/env Rscript

# Create the directory for analyses
if (!dir.exists(file.path("data","analyses"))) dir.create(file.path("data","analyses"))

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(tidygraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(fluxweb))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(abind))

# Check if bootstrap_forage_ratio exist, otherwise generate it.
if(!file.exists(file.path("data", "processed", "bootstrap_forage_ratio.csv"))) {
  system(paste("nohup Rscript", file.path("code", "ModelForageResponse.R")))
  }

# Source all needed functions for the analyses
source("./code/CalculateFluxes.R")

# Import the data
temperature <- read_csv(file = file.path("data", "processed", "interpolation", "temperature.csv"), show_col_types = FALSE)
weekly_biomasses <- read_csv(file = file.path("data", "processed", "interpolation", "weekly_biomasses.csv"), show_col_types = FALSE)
weekly_bodymass <- read_csv(file = file.path("data", "processed", "interpolation", "weekly_bodymass.csv"), show_col_types = FALSE)
node_data <- read_csv(file = file.path("data", "raw", "node_data.csv"), show_col_types = FALSE)
forage_ratio <- read_csv(file = file.path("data", "processed", "forage_ratio.csv"), show_col_types = FALSE)
bootstrap_forage_ratio <- read_csv(file = file.path("data", "processed", "bootstrap_forage_ratio.csv"), show_col_types = FALSE)
  
# Define the stations and the dates we want to compute the fluxes
station = "BY31 LANDSORTSDJ"
dates <- weekly_biomasses |>
  # Filter weekly_biomasses to only include samples from 2007 to end of 2023
  filter(sample_week > as.Date("2007-03-05"),
         sample_week <= as.Date("2023-12-04")) |> 
  pull(sample_week) |> unique()

# Null model ----
# 
# Run the presence_absence model and save it as the null model
if(!file.exists(file.path("data", "analyses", "timeseries_null.csv"))){
  plan(multisession)
  presence_absence_model <- tibble(sample_week = dates) |> 
    mutate(flux_df = future_map(sample_week, ~ {
      flux_mat <- possibly(dna2flux, otherwise = NULL)(
        forage_ratio = forage_ratio,
        node_data = node_data,
        weekly_biomasses = weekly_biomasses,
        weekly_bodymass = weekly_bodymass,
        temperature = temperature,
        date = .x,
        station = station,
        as_graph = FALSE,
        presence_absence = TRUE
      )
      if(!is.null(flux_mat)) {
        flux_mat |> 
          as.data.frame() |> 
          rownames_to_column(var = "predator") |> 
          pivot_longer(-predator, names_to = "prey", values_to = "flux") |> 
          filter(flux > 0)
      }
    }, .options = furrr_options(seed = TRUE))) |> 
    unnest(flux_df) |> 
    mutate(station = station)
  write_csv(presence_absence_model, file.path("data", "analyses", "timeseries_null.csv"))
}

# Timeseries model ----
#
# Calculate fluxes using the forage ratio based on the 1000 iterations for each date...
# But it takes some time, so it is better to do the computations once and save for later if we need to reuse

# Define the directory we want to save the bootrapped fluxes to be save
cache.dir = file.path("data", "analyses", "bootstrapped_fluxes")

# Parallelize for faster computations
plan(multisession)
future_walk(dates, function(date) {
  cacheMyFluxes(cache.dir = cache.dir,
                bootstrap_forage_ratio = bootstrap_forage_ratio,
                node_data = node_data,
                weekly_biomasses = weekly_biomasses,
                weekly_bodymass = weekly_bodymass,
                temperature = temperature,
                date = date,
                station = station,
                as_graph = FALSE)
  })

# Create a tibble containing all dates and graph object with confidence interval and save it for later
if(!file.exists(file.path("data", "analyses", "as_tbl_graph_timeseries_fluxes.rds"))){
  plan(multisession)
  timeseries_fluxes <- tibble(sample_week = dates) |>
    # For each unique sampling week, compute the confidence graph in parallel
    mutate(conf_graph = future_map(sample_week, function(sample_date) {
      fluxingWithConfidence(bootstrap_forage_ratio = bootstrap_forage_ratio,
                            node_data = node_data,
                            weekly_biomasses = weekly_biomasses, 
                            weekly_bodymass = weekly_bodymass,
                            temperature = temperature,
                            date = sample_date,
                            station = station,
                            cache.dir = cache.dir) |> 
        # Switch to node view (for adding metadata)
        activate(nodes) |>
        # Add sample week
        mutate(sample_week = sample_date)
    }, .options = furrr_options(seed = TRUE)  # Ensures reproducibility
    ))
  # Save the daily_fluxes tibble as a rds file
  write_rds(timeseries_fluxes, file = file.path("data", "analyses", "as_tbl_graph_timeseries_fluxes.rds"))
}

# Modify the tibble to save timeseries_fluxes as a dataframe with sample_week, year, iso_week, predator, prey, mean +- CI fluxes
if(!file.exists(file.path("data", "analyses", "timeseries_fluxes.csv"))){
  plan(multisession)
  timeseries_df <- timeseries_fluxes |>
    pull(conf_graph) |> 
    future_map_dfr(~ .x |>
                     activate(edges) |>
                     mutate(across(c(flux_mean, flux_lower, flux_upper), ~ na_if(.x, 0))) |>
                     extract_flux_long(),
                   .options = furrr_options(seed = TRUE)) |>
    mutate(iso_week = isoweek(sample_week),
           year = year(sample_week),
           station = station) |>
    select(sample_week, year, iso_week,station, predator, prey, mean = flux_mean, upper = flux_upper, lower = flux_lower)
  # Save as csv    
  write_csv(timeseries_df, file = file.path("data", "analyses", "timeseries_fluxes.csv"))
}

# Time aggregated fluxes ----
 
## Aggregate per iso_week ----
if(!file.exists(file.path("data", "analyses", "isoweek_fluxes.csv"))){
  # Create a tibble that have all the paths for each files
  files <- tibble(path = list.files(cache.dir, full.names = TRUE)) |>
    mutate(filename = basename(path),
           date = str_extract(filename, "\\d{4}-\\d{2}-\\d{2}") |> as.Date(),
           year = year(date),
           iso_week = isoweek(date),
           station = str_extract(filename, "(?<=flux_)[^_]+"))
  
  plan(multisession)
  isoweek_fluxes <- files |>
    group_nest(iso_week, station) |> 
    mutate(flux = future_map(data, ~ {
      aggregated_fluxes <- aggregateFluxes(.x$path)
      if (is.null(aggregated_fluxes)) return(NULL)
      aggregated_fluxes
    }, .options = furrr_options(seed = TRUE)
    )) |> 
    select(- data) |> 
    unnest(flux) |> 
    filter(mean > 0 & lower > 0 & upper > 0)
  # Save as csv  
  write_csv(isoweek_fluxes, file = file.path("data", "analyses", "isoweek_fluxes.csv"))
}

## Aggregate per year ----
if(!file.exists(file.path("data", "analyses", "annual_fluxes.csv"))){
  files <- tibble(path = list.files(cache.dir, full.names = TRUE)) |>
    mutate(filename = basename(path),
           date = str_extract(filename, "\\d{4}-\\d{2}-\\d{2}") |> as.Date(),
           year = year(date),
           iso_week = isoweek(date),
           station = str_extract(filename, "(?<=flux_)[^_]+"))
  plan(multisession)
  annual_fluxes <- files |>
    filter(year > 2007) |> # As 2007 starts with values in March and misses January and February
    group_nest(year, station) |>
    mutate(flux = future_map(data, ~ {
      aggregated_fluxes <- aggregateFluxes(.x$path)
      if (is.null(aggregated_fluxes)) return(NULL)
      aggregated_fluxes
    }, .options = furrr_options(seed = TRUE)
    )) |>
    select(- data) |> 
    unnest(flux) |> 
    filter(mean > 0 & lower > 0 & upper > 0)
  # Save as csv  
  write_csv(annual_fluxes, file = file.path("data", "analyses", "annual_fluxes.csv"))
}

# Visualisations ----
color_mapping = setNames(node_data$color, node_data$node_name)

## GIF ----
if(!file.exists(file.path("output", "figure", "dna_flux.gif"))){
  # Plotting helper function
  plot_fluxes <- function(data, trophic_level, unit = "kJ/day/m2", abline = FALSE) {
    abline_x <- if (!isFALSE(abline)) as.numeric(as.Date(abline) - as.Date("1970-01-01")) else NULL
    
    p <- data |>
      filter(trophic_level == !!trophic_level) |>
      ggplot(aes(x = sample_date, y = flux_mean, ymin = flux_lower, ymax = flux_upper)) +
      geom_line(aes(color = prey)) +
      geom_ribbon(aes(fill = prey), alpha = 0.4) +
      facet_grid(predator ~ .) +
      scale_fill_manual(values = color_mapping) +
      scale_color_manual(values = color_mapping) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(y = paste0("Fluxes (", unit, ")"), x = NULL)
    
    if (!is.null(abline_x)) {
      p <- p + geom_vline(xintercept = abline_x, linetype = "dashed", color = "red")
    }
    
    return(p)
  }
  
  # Processing data helper function
  fluxTimeSeries <- function(file, time_unit, multiplier = 1, unit_label, abline = FALSE) {
    data <- read_csv(file.path("data", "analyses", file), show_col_types = FALSE) |>
      mutate(flux_mean = mean * multiplier,
             flux_lower = lower * multiplier,
             flux_upper = upper * multiplier,
             sample_date = sample_week) |> 
      left_join(node_data |> rename(predator = node_name), by = "predator") |>
      mutate(prey = factor(prey, levels = node_data$node_name))
    
    cowplot::plot_grid(
      plot_fluxes(data, trophic_level = 2, unit = unit_label, abline = abline),
      plot_fluxes(data, trophic_level = 3, unit = unit_label, abline = abline) +
        theme(axis.text.x = element_text(), axis.ticks.x = element_line()),
      ncol = 1,
      rel_heights = c(8, 3),
      align = "v"
    )
  }
  
  dir.create(file.path("output", "figure", "gif_frames"), showWarnings = FALSE)
  df <- read_rds(file.path("data", "analyses", "as_tbl_graph_timeseries_fluxes.rds")) |> 
    pull(conf_graph) |>
    map(~ .x |>
          activate(edges) |> 
          mutate(flux_mean = ifelse(flux_mean == 0, NA, flux_mean),
                 flux_lower = ifelse(flux_lower == 0, NA, flux_lower),
                 flux_upper = ifelse(flux_upper == 0, NA, flux_upper)))
  

  edge_color_gradient <- scale_edge_color_gradientn(
    colours = c("grey50", "#A23C2A", "#136F63"),
    values = scales::rescale(c(log10(1), log10(1.5), log10(5)), from = c(log10(1), log10(5))),
    limits = c(log10(1), log10(5))
  )
  edge_width_scale <- scale_edge_width(range = c(1, 20), limits = c(log10(1), log10(5)))
  fill_scale <- scale_fill_manual(values = color_mapping)
  theme_options <- theme_graph() + theme(
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_blank()
  )
    

  plan(multisession)
  # Use future_map to process frames in parallel
  future_map(1:length(df), function(i) {

    week_data <- df[[i]]
    
    week <- week_data |>
      activate(nodes) |>
      as_tibble() |>
      pull(sample_week) |>
      unique()
    
    ts <- fluxTimeSeries("timeseries_fluxes.csv",
                         multiplier = 1,
                         unit_label = "kJ/day/m2",
                         abline = as.Date(week))
    # Create the plot
    foodweb <- week_data |>
      ggraph(layout = "manual", x = horizontal_position, y = trophic_level) +
      geom_edge_link(aes(edge_width = log10(flux_mean + 1), col = log10(flux_mean + 1)),
                     arrow = arrow(length = unit(2, 'mm'), ends = "first")) +
      geom_node_point(aes(size = pi*(biomass)^2, fill = name), shape = 21) +
      geom_node_label(aes(label = name), angle = 45, size = 3, nudge_y = -0.05, nudge_x = .2, hjust = 1) +
      edge_color_gradient +
      edge_width_scale +
      fill_scale +
      theme_options +
      coord_cartesian(xlim = c(0,14),
                      ylim = c(0,3.5)) +
      annotate("label", x = 2.8, y = 3.2,
               label = paste("Week:", week), size = 5,
               label.size = 0.5, label.r = unit(0.15, "lines"),
               fill = "#ffffffcc", color = "black", fontface = "bold")
    
    plot <- cowplot::plot_grid(ts, foodweb)

    # Save the plot
    ggsave(sprintf("./output/figure/gif_frames/frame_%03d.png", i), plot = plot, width = 16, height = 9, dpi = 300)
  },
  .options = furrr_options(seed = TRUE)
  )

  png_files <- list.files("./output/figure/gif_frames", pattern = "*.png", full.names = TRUE)
  gifski::gifski(png_files, gif_file = "./output/figure/dna_flux.gif", width = 4800, height = 2700, delay = 0.1)
  
  
  unlink(file.path("output", "figure", "gif_frames"), recursive = TRUE)
}

## Impact of forage ratios ----
# or how the forage ratios deviate from the null model
timeseries_fluxes <- read_csv(file.path("data", "analyses", "timeseries_fluxes.csv"), show_col_types = FALSE)
isoweek_fluxes <- read_csv(file.path("data", "analyses", "isoweek_fluxes.csv"), show_col_types = FALSE)
timeseries_null <- read_csv(file.path("data", "analyses", "timeseries_null.csv"), show_col_types = FALSE)
isoweek_null <- timeseries_null |> 
  group_by(predator, prey, iso_week = isoweek(sample_week), station) |> 
  summarise(avg_fluxes = mean(flux, na.rm = T), .groups = "drop")

timeseries_diff <-
  timeseries_fluxes |> 
  left_join(timeseries_null |> rename("null" = flux),
            by = c("sample_week","station", "predator", "prey")) |> 
  mutate(predator = factor(predator, levels = node_data$node_name),
         prey = factor(prey, levels = node_data$node_name)) |> na.omit() |> 
  ggplot(mapping = aes(x = sample_week, y = mean-null, ymin = lower-null, ymax = upper-null, col = prey, fill = prey))+
  geom_hline(yintercept = 0)+
  geom_ribbon(alpha = .2,linewidth = .2) +
  geom_line()+
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  facet_grid(predator~., scales = "free") +
  scale_x_date(date_breaks = "4 months", expand = c(0, 0))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y = "Deviance from model without selectivity \n [kJ/day/m2]",
       x = NULL)
ggsave(plot = timeseries_diff,
       filename = file.path("output", "figure", "timeseries_diff.pdf"),
       width = 12,
       height = 10)

# Reusable plot function
plot_flux_difference <- function(trophic_level_filter) {
  isoweek_fluxes |> 
    left_join(isoweek_null |> rename("null" = avg_fluxes),
              by = c("iso_week", "station", "predator", "prey")) |> 
    mutate(
      predator = factor(predator, levels = node_data$node_name),
      prey = factor(prey, levels = node_data$node_name),
      trophic_level = case_when(
        prey %in% node_data$node_name[node_data$trophic_level == 2] ~ "zooplankton",
        prey %in% node_data$node_name[node_data$trophic_level == 1] ~ "phytoplankton"
      )
    ) |>
    filter(trophic_level == trophic_level_filter,
           iso_week %in% 2:52) |>
    ggplot(aes(x = iso_week)) +
    
    geom_ribbon(aes(y = mean, ymin = lower, ymax = upper),
                alpha = .2, linewidth = .2,
                fill = "black", col = "black") +
    
    geom_line(aes(y = mean, color = "Selectivity")) +
    geom_line(aes(y = null, color = "No selectivity")) +
    
    facet_grid(prey ~ predator, scales = "free_y") +
    
    scale_color_manual(
      name = NULL,
      values = c("Selectivity" = "black", "No selectivity" = "#ff7f00")
    ) +
    
    scale_x_continuous(
      breaks = seq(1, 52.1775, 4.348125),
      labels = month.abb, expand = c(0, 0)
    ) +
    
    theme_bw() +
    labs(
      y = "Fluxes \n [kJ/day/m2]",
      x = NULL
    )
}

# Generate the two plots
(phytoplankton_diff <- plot_flux_difference("phytoplankton"))
(zooplankton_diff   <- plot_flux_difference("zooplankton"))


isoweek_rel_diff <-
  isoweek_fluxes |> 
  left_join(isoweek_null |> rename("null" = avg_fluxes),
            by = c("iso_week", "station", "predator", "prey")) |> 
  group_by(iso_week, predator) |> 
  mutate(predator = factor(predator, levels = node_data$node_name),
         prey = factor(prey, levels = node_data$node_name),
         delta_mean = mean / sum(mean, na.rm = TRUE) * 100,
         delta_upper = upper / sum(mean, na.rm = TRUE) * 100,
         delta_lower = lower / sum(mean, na.rm = TRUE) * 100,
         delta_null = null / sum(null, na.rm = TRUE) * 100) |> na.omit() |> 
  ungroup() |> 
  
  ggplot(mapping = aes(x = iso_week, y = delta_mean-delta_null, ymin = delta_lower-delta_null, ymax = delta_upper-delta_null, col = prey, fill = prey))+
  geom_hline(yintercept = 0)+
  geom_ribbon(alpha = .2,linewidth = .2) +
  geom_line()+
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  facet_grid(predator~., scales = "free") +
  scale_x_continuous(breaks = seq(1,52.1775, 4.348125),
                     labels = month.abb, expand = c(0, 0))+
  theme_bw() +
  labs(y = "Deviance from model without selectivity \n [%]",
       x = NULL)
ggsave(plot = isoweek_rel_diff,
       filename = file.path("output", "figure", "isoweek_rel_diff.pdf"),
       width = 6,
       height = 12)

## Predation pressure ----

# Load data
weekly_biomasses <- read_csv("data/processed/interpolation/weekly_biomasses.csv", show_col_types = FALSE)
timeseries_fluxes <- read_csv("data/analyses/timeseries_fluxes.csv", show_col_types = FALSE)
isoweek_fluxes <- read_csv("data/analyses/isoweek_fluxes.csv", show_col_types = FALSE)

# Join and calculate predation pressure
timeseries_pressure <- weekly_biomasses |> 
  filter(
    station_name == "BY31 LANDSORTSDJ",
    !node_name %in% c("Clupea", "Gasterosteus", "Sprattus")
  ) |> 
  rename_with(~c("prey", "station"), .cols = c("node_name", "station_name")) |> 
  select(prey, station, sample_week, biomass) |> 
  right_join(timeseries_fluxes, by = c("prey", "station", "sample_week"), relationship = "many-to-many") |> 
  mutate(
    iso_week = isoweek(sample_week),
    year = as.factor(year(sample_week))
  ) |>
  group_by(prey, sample_week, station, iso_week, year, biomass) |> 
  summarise(pressure_sum = sum(mean, na.rm = TRUE), .groups = "drop_last") |> 
  summarise(pressure = sum(pressure_sum, na.rm = TRUE) / sum(biomass, na.rm = TRUE), .groups = "drop") |> 
  filter(iso_week %in% 2:52) |> 
  mutate(
    trophic_level = case_when(
      prey %in% node_data$node_name[node_data$trophic_level == 2] ~ "zooplankton",
      prey %in% node_data$node_name[node_data$trophic_level == 1] ~ "phytoplankton"
    )
  )

# Add total predation pressure across prey
timeseries_pressure_with_tot <- timeseries_pressure |> 
  group_by(prey = "Total", sample_week, station, iso_week, year, trophic_level) |> 
  summarise(pressure = sum(pressure, na.rm = TRUE), .groups = "drop") |> 
  bind_rows(timeseries_pressure) |> 
  mutate(prey = factor(prey, levels = c("Total", node_data$node_name)))

# Aggregate by isoweek
isoweek_pressure <- timeseries_pressure_with_tot |> 
  group_by(prey, iso_week, station, trophic_level) |> 
  summarise(pressure = median(pressure, na.rm = TRUE), .groups = "drop")

# Plotting function
plot_pressure <- function(tl, label_y = TRUE) {
  ggplot(mapping = aes(x = iso_week, y = pressure)) +
    geom_point(
      data = filter(timeseries_pressure_with_tot, trophic_level == tl),
      mapping = aes(group = year), size = 0.2
    ) +
    geom_line(
      data = filter(isoweek_pressure, trophic_level == tl),
      col = "#C44536", linewidth = 0.8
    ) +
    facet_grid(prey ~ ., scales = "fixed") +
    theme_bw() +
    scale_x_continuous(
      breaks = seq(1, 52.1775, 4.348125),
      labels = month.abb, expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = if (label_y) "Predation pressure \n [kJ/day/g]" else NULL
    )
}

# Generate plots
pressure_phyto <- plot_pressure("phytoplankton", label_y = TRUE)
pressure_zoo <- plot_pressure("zooplankton", label_y = FALSE)
blank_spacer <- ggplot() + theme_void()
# Combine
library(patchwork)
predation_pressure_plot <-
  pressure_phyto + 
  (pressure_zoo / blank_spacer + plot_layout(heights = c(8, 13-8))) + 
  plot_layout(widths = c(1, 1))
ggsave(plot = predation_pressure_plot, filename = file.path("output", "figure", "predation_pressure.pdf"), width =10, height = 12)


# Join and calculate predation pressure
timeseries_pred_pressure <- weekly_biomasses |> 
  filter(
    station_name == "BY31 LANDSORTSDJ",
    !node_name %in% c("Clupea", "Gasterosteus", "Sprattus")
  ) |> 
  rename_with(~c("prey", "station"), .cols = c("node_name", "station_name")) |> 
  select(prey, station, sample_week, biomass) |> 
  right_join(timeseries_fluxes, by = c("prey", "station", "sample_week"), relationship = "many-to-many") |> 
  mutate(
    iso_week = isoweek(sample_week),
    year = as.factor(year(sample_week))
  ) |>
  group_by(prey, predator, sample_week, station, iso_week, year) |> 
  summarise(pressure = sum(mean, na.rm = TRUE) / sum(biomass, na.rm = TRUE), .groups = "drop") |> 
  filter(iso_week %in% 2:52) |> 
  mutate(
    trophic_level = case_when(
      prey %in% node_data$node_name[node_data$trophic_level == 2] ~ "zooplankton",
      prey %in% node_data$node_name[node_data$trophic_level == 1] ~ "phytoplankton"
    )
  )

# Add total predation pressure across prey
timeseries_pred_pressure_with_tot <- timeseries_pred_pressure |> 
  group_by(prey = "Total",predator, sample_week, station, iso_week, year, trophic_level) |> 
  summarise(pressure = sum(pressure, na.rm = TRUE), .groups = "drop") |> 
  bind_rows(timeseries_pred_pressure) |> 
  mutate(prey = factor(prey, levels = c("Total", node_data$node_name)),
         predator = factor(predator, levels = c(node_data$node_name))) |> filter(prey == "Total")

# Aggregate by isoweek
isoweek_pred_pressure <- timeseries_pred_pressure_with_tot |> 
  group_by(prey, predator, iso_week, station, trophic_level) |> 
  summarise(median_pressure = median(pressure, na.rm = TRUE),
            pressure_upper = quantile(pressure, 0.75, na.rm = TRUE),
            pressure_lower = quantile(pressure, 0.25, na.rm = TRUE), .groups = "drop")

ggplot(data = isoweek_pred_pressure,mapping = aes(x = iso_week, y = median_pressure, ymin = pressure_lower, ymax = pressure_upper,  col = predator, fill = predator)) +
  geom_ribbon(alpha = .1)+
  geom_line(linewidth = 1) +
  facet_grid(trophic_level~., scales = "free") +
  theme_bw() +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb, expand = c(0, 0)
  ) +
  scale_color_manual(values = color_mapping)+
  scale_fill_manual(values = color_mapping)+
  
    labs(
    x = NULL,
    y = "Predation pressure \n (median +- interquartile ranges) \n [kJ/day/g]"
  )
