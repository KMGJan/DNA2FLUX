#!/usr/bin/env Rscript

## Load data
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(tidygraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(rlang))
source("./code/CalculateFluxes.R")

if (!dir.exists(file.path("data","analyses"))) dir.create(file.path("data","analyses"))

if(!file.exists(file.path("data", "analyses", "weekly_fluxes.rds"))){
  # Load the needed libraries
  suppressPackageStartupMessages(library(fluxweb))
  suppressPackageStartupMessages(library(furrr))
  
  # Import the data
  temperature <- read_csv(file = file.path("data", "processed", "interpolation", "temperature.csv"), show_col_types = FALSE)
  weekly_biomasses <- read_csv(file = file.path("data", "processed", "interpolation", "weekly_biomasses.csv"), show_col_types = FALSE)
  weekly_bodymass <- read_csv(file = file.path("data", "processed", "interpolation", "weekly_bodymass.csv"), show_col_types = FALSE)
  node_data <- read_csv(file = file.path("data", "raw", "node_data.csv"), show_col_types = FALSE)
  forage_ratio <- read_csv(file = file.path("data", "processed", "forage_ratio.csv"), show_col_types = FALSE)
  bootstrap_forage_ratio <- read_csv(file = file.path("data", "processed", "bootstrap_forage_ratio.csv"), show_col_types = FALSE)
  
  # Define the stations and the dates we want to compute
  station = "BY31 LANDSORTSDJ"
  dates <-
    weekly_biomasses |>
    # Filter weekly_biomasses to only include samples from 2007 to end of 2023
    filter(sample_week > as.Date("2007-03-05"),
           sample_week <= as.Date("2011-12-04")) |> 
    pull(sample_week) |> unique()

  
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
  # Do the weekly calculations in parallel
  plan(multisession)
  tibble(sample_week = dates) |>
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
        )) #|> 
      # Save the final tibble as a rds file
      write_rds(file = file.path("data", "analyses", "weekly_model.rds"))
 }


file_tbl <- list.files(cache.dir, pattern = glue::glue("flux_{station}_.*\\.rds$"), full.names = TRUE) |>
  tibble(path = _) |>
  mutate(
    filename = basename(path),
    date = str_extract(filename, "\\d{4}-\\d{2}-\\d{2}") |> as.Date(),
    month = floor_date(date, "month"),
    year = year(date)
  )

# --- Step 2: Aggregate by month ---
plan(multisession)  # or multicore, depending on your OS
library(abind)
monthly_fluxes <- file_tbl |>
  group_nest(month, year) |>
  mutate(
    flux = future_map(data, ~ aggregateFluxes(.x$path), 
                      .options = furrr_options(seed = TRUE))
  ) |>
  select(year, month, flux) |>
  filter(!map_lgl(flux, is.null))

# Optional: unnest the result to get one row per predator-prey-month
monthly_fluxes_long <- monthly_fluxes |> 
  unnest(flux) |> 
  filter(mean > 0 & lower > 0 & upper > 0)


monthly_fluxes_long |> 
  mutate(mean = (mean) * 365,
         lower = lower * 365,
         upper = upper * 365) |> 
  left_join(node_data |>
              rename("predator" = node_name),
            by = "predator") |> 
  filter(trophic_level == 3) |> 
  ggplot(aes(x = month, y = mean+1, ymin = lower+1, ymax = upper+1)) +
  geom_line(mapping = aes(col = prey)) +
  geom_ribbon(mapping = aes(fill = prey), alpha = .4)+
  facet_grid(predator~.)+
  #scale_fill_manual(values = color_mapping)+
  #scale_color_manual(values = color_mapping) +
  theme_bw()+
  scale_y_log10()+
  annotation_logticks() +
  labs(y = "Fluxes (kJ/d/m2)",
       x = NULL)

rds_files <- list.files(cache.dir, pattern = "\\.rds$", full.names = TRUE)

# Read each .rds file and extract its dimensions
dim_report <- rds_files %>%
  map(~ {
    dims <- dim(read_rds(.))
    tibble(file = basename(.), 
           dim_1 = dims[1], 
           dim_2 = dims[2], 
           dim_3 = dims[3],
           date = str_extract(basename(.), "\\d{4}-\\d{2}-\\d{2}"))  # Extract date in YYYY-MM-DD format
  }) %>%
  bind_rows()

# View the report
print(dim_report)
dim_report |>
 # filter(dim_3 < 250) |> 
  mutate(date = as.Date(date)) |> 
  ggplot(aes(x = date, y = dim_3))+
  geom_line()

weekly_model <- read_rds(file = file.path("data", "processed", "weekly_model.rds"))
color_mapping = setNames(node_data$color, node_data$node_name)
df <- weekly_model |> 
  pull(conf_graph) |>
  map(~ .x |>
        activate(edges) |> 
        mutate(flux_mean = ifelse(flux_mean == 0, NA, flux_mean),
               flux_lower = ifelse(flux_lower == 0, NA, flux_lower),
               flux_upper = ifelse(flux_upper == 0, NA, flux_upper)))


if(!file.exists(file.path("output", "figure", "dna_flux.gif"))){
  suppressPackageStartupMessages(library(furrr))
  dir.create(file.path("output", "figure", "gif_frames"), showWarnings = FALSE)
  df <- weekly_model |> 
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




# Helper function to convert matrix to long format
adj_to_long <- function(graph, attr, value_name) {
  as.matrix(as_adj(graph, attr = attr, sparse = FALSE)) |>
    as.data.frame() |>
    rownames_to_column("predator") |>
    pivot_longer(-predator, names_to = "prey", values_to = value_name)
}

# Wrapper to extract week + long-form fluxes from a tbl_graph
extract_flux_long <- function(graph) {
  week <- graph |>
    activate(nodes) |>
    as_tibble() |>
    pull(sample_week) |>
    unique()
  
  # Get all three matrices as long-form data frames
  adj_long_mean <- adj_to_long(graph, "flux_mean", "flux_mean")
  adj_long_upper <- adj_to_long(graph, "flux_upper", "flux_upper")
  adj_long_lower <- adj_to_long(graph, "flux_lower", "flux_lower")
  
  # Join and annotate
  flux_long <- adj_long_mean |>
    filter(!is.na(flux_mean)) |>
    left_join(adj_long_upper, by = c("predator", "prey")) |>
    left_join(adj_long_lower, by = c("predator", "prey")) |>
    mutate(sample_week = week)
  
  return(flux_long)
}

all_fluxes <- map_dfr(df, extract_flux_long)
all_fluxes$flux_mean |> min()
p1 <- all_fluxes |> 
  left_join(node_data |>
              rename("predator" = node_name),
            by = "predator") |> 
  filter(trophic_level == 2) |> 
  ggplot(aes(x = sample_week, y = flux_mean+1, ymin = flux_lower+1, ymax = flux_upper+1)) +
  geom_line(mapping = aes(col = prey)) +
  geom_ribbon(mapping = aes(fill = prey), alpha = .4)+
  facet_grid(predator~.)+
  scale_fill_manual(values = color_mapping)+
  scale_color_manual(values = color_mapping) +
  theme_bw()+
  scale_y_log10()+
  annotation_logticks() +
  labs(y = "Fluxes (kJ/week/m2)",
       x = NULL)
p2 <- all_fluxes |> 
  left_join(node_data |>
              rename("predator" = node_name),
            by = "predator") |> 
  filter(trophic_level == 3) |> 
  ggplot(aes(x = sample_week, y = flux_mean+1, ymin = flux_lower+1, ymax = flux_upper+1)) +
  geom_line(mapping = aes(col = prey)) +
  geom_ribbon(mapping = aes(fill = prey), alpha = .4)+
  facet_grid(predator~.)+
  scale_fill_manual(values = color_mapping)+
  scale_color_manual(values = color_mapping) +
  theme_bw()+
  scale_y_log10()+
  annotation_logticks() +
  labs(y = "Fluxes (kJ/week/m2)",
       x = NULL)
timeseries<-cowplot::plot_grid(p1,p2, ncol = 1, rel_heights = c(8,3), align = T)
ggsave(plot = timeseries,
       filename = file.path("output", "figure", "weekly_timeseries.pdf"),
       height = 10,
       width = 10)
test <- aggregateFlux(cache.dir = cache.dir, aggregation_period = "yearly")

