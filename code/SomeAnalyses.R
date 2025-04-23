source("./code/CalculateFluxes.R")
## Load data
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(fluxweb))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(tidygraph))
suppressPackageStartupMessages(library(igraph))

temperature <- read_csv(file = file.path("data", "processed", "interpolation", "temperature.csv"), show_col_types = FALSE)
weekly_biomasses <- read_csv(file = file.path("data", "processed", "interpolation", "weekly_biomasses.csv"), show_col_types = FALSE)
weekly_bodymass <- read_csv(file = file.path("data", "processed", "interpolation", "weekly_bodymass.csv"), show_col_types = FALSE)
node_data <- read_csv(file = file.path("data", "raw", "node_data.csv"), show_col_types = FALSE)
forage_ratio <- read_csv(file = file.path("data", "processed", "forage_ratio.csv"), show_col_types = FALSE)
bootstrap_forage_ratio <- read_csv(file = file.path("data", "processed", "bootstrap_forage_ratio.csv"), show_col_types = FALSE) |> #For 200 bootstraps it takes approx. 20 min to run all weeks
  filter(Iteration <= 1000)

# Run
plan(multisession)
dates_tbl <- weekly_biomasses |>
  filter(year(sample_week) >= 2007,
         sample_week <= as.Date("2023-12-04")) |> 
  distinct(sample_week)

conf_tbl <- dates_tbl |>
  mutate(conf_graph = future_map(sample_week, function(j) {
    g <- fluxConfidence(
      bootstrap_forage_ratio, node_data, weekly_biomasses,
      weekly_bodymass, temperature,
      date = j,
      station = "BY31 LANDSORTSDJ",
      parallel = FALSE
    )
    if (!is.null(g)) {
      g <- activate(g, nodes) |> mutate(sample_week = j)
    }
    g
  },
  .options = furrr_options(seed = TRUE)  # Ensures reproducibility
  ))

conf_tbl |> 
  write_rds(file = file.path("data", "processed", "weekly_model.rds"))

weekly_model <- read_rds(file = file.path("data", "processed", "weekly_model.rds"))

df <- weekly_model |> 
  #filter(sample_week == "2007-03-19") |> 
  pull(conf_graph) |> 
  map(~ .x |> activate(edges) |> mutate(mean = ifelse(mean == 0, NA, mean),
                                        lower_ci = ifelse(lower_ci == 0, NA, lower_ci),
                                        upper_ci = ifelse(upper_ci == 0, NA, upper_ci)))


color_mapping = setNames(node_data$color, node_data$node_name)
dir.create("./output/figure/gif_frames", showWarnings = FALSE)
for(i in 1:length(df)){
  week <- df[[i]] |>
    activate(nodes) |>
    as_tibble() |>
    pull(sample_week) |>
    unique()
  
  plot <- df[[i]] |>
    ggraph(layout = "manual", x = horizontal_position, y = trophic_level) +
    geom_edge_link(aes(edge_width = mean, col = mean),
                   arrow = arrow(length = unit(3, 'mm'), ends = "first")) +
    geom_node_point(aes(size = pi*(biomass)^2, fill = name), shape = 21) +
    geom_node_label(aes(label = name), angle = 45, size = 3, nudge_y = -0.05) +
    theme_graph() +
    scale_fill_manual(values = color_mapping)+
    scale_edge_color_gradient2(
      low = "#ffffbf",
      mid = "#2c3d2b",
      high = "#fc8d59",
      midpoint = 0.05) +
    theme(legend.position = "none") +
    ggtitle(paste("Week:", week))
  
  ggsave(sprintf("./output/figure/gif_frames/frame_%03d.png", i), plot = plot, width = 8, height = 10, dpi = 150)
}
png_files <- list.files("./output/figure/gif_frames", pattern = "*.png", full.names = TRUE)
gifski::gifski(png_files, gif_file = "./output/figure/dna_flux.gif", width = 800, height = 1000, delay = 0.3)



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
  adj_long_mean <- adj_to_long(graph, "mean", "flux_mean")
  adj_long_upper <- adj_to_long(graph, "upper_ci", "flux_upper")
  adj_long_lower <- adj_to_long(graph, "lower_ci", "flux_lower")
  
  # Join and annotate
  flux_long <- adj_long_mean |>
    filter(!is.na(flux_mean)) |>
    left_join(adj_long_upper, by = c("predator", "prey")) |>
    left_join(adj_long_lower, by = c("predator", "prey")) |>
    mutate(sample_week = week)
  
  return(flux_long)
}
flux_df <- extract_flux_long(df[[1]])
all_fluxes <- map_dfr(df, extract_flux_long)

