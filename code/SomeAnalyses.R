#!/usr/bin/env Rscript

# Create the directory for analyses
if (!dir.exists(file.path("data","analyses"))) dir.create(file.path("data","analyses"))

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(tidygraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(rlang))

# Source all needed functions for the analyses
source("./code/CalculateFluxes.R")

# If the daily_fluxes.rds tibble was not produced before, this will run
if(!file.exists(file.path("data", "analyses", "as_tbl_graph_daily_fluxes.rds")) |
   !file.exists(file.path("data", "analyses", "daily_fluxes.csv")) |
   !file.exists(file.path("data", "analyses", "monthly_fluxes.csv")) |
   !file.exists(file.path("data", "analyses", "yearly_fluxes.csv"))){
  
  # Load the needed libraries
  suppressPackageStartupMessages(library(fluxweb))
  suppressPackageStartupMessages(library(furrr))
  suppressPackageStartupMessages(library(abind))
  
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
    pull(sample_week) |>
    unique()
  
  # As bootstrapping with 1000 iterations is time consuming (ca. 2h)
  # It is better to do the computations once and save for later if we need to reuse
  
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

  # Create a tibble containing all dates and graph object with confidence interval 
  plan(multisession)
  daily_fluxes <- tibble(sample_week = dates) |>
    # For each unique sampling week, compute the confidence graph in parallel
    mutate(conf_graph = future_map(sample_week, function(sample_date) {
      fluxingWithConfidence(
          bootstrap_forage_ratio = bootstrap_forage_ratio,
          node_data = node_data,
          weekly_biomasses = weekly_biomasses, 
          weekly_bodymass = weekly_bodymass,
          temperature = temperature,
          date = sample_date,
          station = station,
          cache.dir = cache.dir
          ) |> 
          
        # Switch to node view (for adding metadata)
        activate(nodes) |>
        # Add sample week
        mutate(sample_week = sample_date)
        },
        .options = furrr_options(seed = TRUE)  # Ensures reproducibility
        ))
  
  # Save the daily_fluxes tibble as a rds file
  write_rds(daily_fluxes, file = file.path("data", "analyses", "as_tbl_graph_daily_fluxes.rds"))

  # Save daily_fluxes as a csv file too
  daily_df <- daily_fluxes |> 
    pull(conf_graph) |> 
    map_dfr(~ .x |> 
              activate(edges) |> 
              mutate(across(c(flux_mean, flux_lower, flux_upper), ~ na_if(.x, 0))) |> 
              extract_flux_long()) |> 
    mutate(iso_week = isoweek(sample_week),
           year = year(sample_week)) |> 
    select(sample_week, year, iso_week, predator, prey, flux_mean, flux_upper, flux_lower)
  
  write_csv(daily_df, file = file.path("data", "analyses", "daily_fluxes.csv"))
    
    # Define a small tibble that contains the month abbreviations and the numeric months
    month_abb <- tibble(month.abb = month.abb,
                        month = 1:12)
    # Create a tibble that have all the paths for each
    files <- tibble(path = list.files(cache.dir, full.names = TRUE)) |>
      mutate(filename = basename(path),
             date = str_extract(filename, "\\d{4}-\\d{2}-\\d{2}") |> as.Date(),
             month = month(date),
             year = year(date),
             station = str_extract(filename, "(?<=flux_)[^_]+")) |> 
      left_join(month_abb, by = "month")
    
    # Aggregate per month
    plan(multisession)
    monthly_fluxes <- files |>
      group_nest(year, month, month.abb, station) |>
      mutate(flux = future_map(
        data, ~ {
          aggregated_fluxes <- aggregateFluxes(.x$path)
          if (is.null(aggregated_fluxes)) return(NULL)
          aggregated_fluxes
          }, .options = furrr_options(seed = TRUE)
        )) |>
      select(- data) |> 
      unnest(flux) |> 
      filter(mean > 0 & lower > 0 & upper > 0)
    
    write_csv(monthly_fluxes, file = file.path("data", "analyses", "monthly_fluxes.csv"))
  
  
    # Aggregate per year
    plan(multisession)
    yearly_fluxes <- files |>
      group_nest(year, station) |>
      mutate(flux = future_map(
        data, ~ {
          aggregated_fluxes <- aggregateFluxes(.x$path)
          if (is.null(aggregated_fluxes)) return(NULL)
          aggregated_fluxes
          }, .options = furrr_options(seed = TRUE)
        )) |>
      select(- data) |> 
      unnest(flux) |> 
      filter(mean > 0 & lower > 0 & upper > 0)
    
    write_csv(yearly_fluxes, file = file.path("data", "analyses", "yearly_fluxes.csv"))
    }

node_data <- read_csv(file = file.path("data", "raw", "node_data.csv"), show_col_types = FALSE)
color_mapping = setNames(node_data$color, node_data$node_name)

monthly_fluxes <- read_csv(file = file.path("data", "analyses", "monthly_fluxes.csv"))
monthly_fluxes |>
  mutate(mean = mean * 30.5,
         lower = lower * 30.5,
         upper = upper * 30.5,
         sample_month = as.Date(paste(year, month, "01", sep = "-"), format = "%Y-%m-%d")) |> 
  left_join(node_data |>
              rename("predator" = node_name),
            by = "predator") |> 
  filter(trophic_level == 2) |>
  ggplot(aes(x = sample_month, y = mean+1, ymin = lower+1, ymax = upper+1)) +
  geom_line(mapping = aes(col = prey)) +
  geom_ribbon(mapping = aes(fill = prey), alpha = .4)+
  facet_grid(predator~.)+
  scale_fill_manual(values = color_mapping)+
  scale_color_manual(values = color_mapping) +
  theme_bw()+
  scale_y_log10()+
  annotation_logticks() +
  labs(y = "Fluxes (kJ/month/m2)",
       x = NULL)
monthly_fluxes |>
  mutate(mean = mean * 30.5,
         lower = lower * 30.5,
         upper = upper * 30.5,
         sample_month = as.Date(paste(year, month, "01", sep = "-"), format = "%Y-%m-%d")) |> 
  left_join(node_data |>
              rename("predator" = node_name),
            by = "predator") |> 
  filter(trophic_level == 3) |>
  ggplot(aes(x = sample_month, y = mean+1, ymin = lower+1, ymax = upper+1)) +
  geom_line(mapping = aes(col = prey)) +
  geom_ribbon(mapping = aes(fill = prey), alpha = .4)+
  facet_grid(predator~.)+
  scale_fill_manual(values = color_mapping)+
  scale_color_manual(values = color_mapping) +
  theme_bw()+
  scale_y_log10()+
  annotation_logticks() +
  labs(y = "Fluxes (kJ/month/m2)",
       x = NULL)

yearly_fluxes <- read_csv(file = file.path("data", "analyses", "yearly_fluxes.csv"))
yearly_fluxes |>
  mutate(mean = mean * 365,
         lower = lower * 365,
         upper = upper * 365) |> 
  left_join(node_data |>
              rename("predator" = node_name),
            by = "predator") |> 
  filter(trophic_level == 2) |>
  ggplot(aes(x = year, y = mean+1, ymin = lower+1, ymax = upper+1)) +
  geom_line(mapping = aes(col = prey)) +
  geom_ribbon(mapping = aes(fill = prey), alpha = .4)+
  facet_grid(predator~.)+
  scale_fill_manual(values = color_mapping)+
  scale_color_manual(values = color_mapping) +
  theme_bw()+
  scale_y_log10()+
  annotation_logticks() +
  labs(y = "Fluxes (kJ/yr/m2)",
       x = NULL)
yearly_fluxes |>
  mutate(mean = mean * 365,
         lower = lower * 365,
         upper = upper * 365) |> 
  left_join(node_data |>
              rename("predator" = node_name),
            by = "predator") |> 
  filter(trophic_level == 3) |>
  ggplot(aes(x = year, y = mean+1, ymin = lower+1, ymax = upper+1)) +
  geom_line(mapping = aes(col = prey)) +
  geom_ribbon(mapping = aes(fill = prey), alpha = .4)+
  facet_grid(predator~.)+
  scale_fill_manual(values = color_mapping)+
  scale_color_manual(values = color_mapping) +
  theme_bw()+
  scale_y_log10()+
  annotation_logticks() +
  labs(y = "Fluxes (kJ/yr/m2)",
       x = NULL)

# Little gif...

if(!file.exists(file.path("output", "figure", "dna_flux.gif"))){
  suppressPackageStartupMessages(library(furrr))
  dir.create(file.path("output", "figure", "gif_frames"), showWarnings = FALSE)
  df <- read_rds(file.path("data", "analyses", "as_tbl_graph_daily_fluxes.rds")) |> 
    pull(conf_graph) |>
    map(~ .x |>
          activate(edges) |> 
          mutate(flux_mean = ifelse(flux_mean == 0, NA, flux_mean),
                 flux_lower = ifelse(flux_lower == 0, NA, flux_lower),
                 flux_upper = ifelse(flux_upper == 0, NA, flux_upper)))
  

  edge_color_gradient <- scale_edge_color_gradientn(
    colours = c("#3D3D3D", "#A23C2A"),
    values = scales::rescale(c(log10(1), log10(5)), from = c(log10(1), log10(5))),
    limits = c(log10(1), log10(5))
  )
  edge_width_scale <- scale_edge_width(range = c(.5, 10), limits = c(log10(1), log10(5)))
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
    # Access df[[i]] just once
    week_data <- df[[i]]
    
    # Extract the unique week once
    week <- week_data |>
      activate(nodes) |>
      as_tibble() |>
      pull(sample_week) |>
      unique()
    
    # Create the plot
    plot <- week_data |>
      ggraph(layout = "manual", x = horizontal_position, y = trophic_level) +
      geom_edge_link(aes(edge_width = log10(flux_mean + 1), col = log10(flux_mean + 1)),
                     arrow = arrow(length = unit(3, 'mm'), ends = "first")) +
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

    # Save the plot
    ggsave(sprintf("./output/figure/gif_frames/frame_%03d.png", i), plot = plot, width = 8.5, height = 7, dpi = 300)
  },
  .options = furrr_options(seed = TRUE)
  )

  png_files <- list.files("./output/figure/gif_frames", pattern = "*.png", full.names = TRUE)
  gifski::gifski(png_files, gif_file = "./output/figure/dna_flux.gif", width = 2550, height = 2100, delay = 0.25)
  
  
  unlink(file.path("output", "figure", "gif_frames"), recursive = TRUE)
}







