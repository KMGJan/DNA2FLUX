#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidygraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(fluxweb))
suppressPackageStartupMessages(library(ggraph))

# Select Station and Day ----

sel_sample_week <- date("2010-07-12")
sel_station_name <- "BY31 LANDSORTSDJ"

read_csv(file = file.path("data", "Processed", "interpolation", "temperature.csv")) |> 
  pull(sample_week) |>  unique()

# Calculate temperature constant
temp <-
  read_csv(file = file.path("data", "Processed", "interpolation", "temperature.csv")) |> 
  filter(sample_week == sel_sample_week,
         station_name == sel_station_name) |> 
  pull(temperature)
boltz <- 0.00008617343  # Boltzmann constant
tkonst <- 0.69/(boltz*(273.15+temp)) #Temperature metabolic constant


# Subset masses for day and station
biomasses <-
  read_csv(file = file.path("data", "processed", "interpolation", "weekly_biomasses.csv")) |> 
  filter(sample_week == sel_sample_week,
         station_name == sel_station_name) |> 
  select(node_name, biomass)

bodymasses <-
  read_csv(file = file.path("data", "processed", "interpolation", "weekly_bodymass.csv")) |> 
  filter(sample_week == sel_sample_week,
         station_name == sel_station_name) |> 
  select(node_name, bodymass)

node_data <-
  read_csv(file = file.path("data", "raw", "node_data.csv")) |> 
  select(node_name, efficiencies, intercept, slope)

# Calculate Weights ----
# (Without bootstrap variation here)

graph <-
  read_csv(file = file.path("data", "processed", "forage_ratio.csv")) |>
  filter(is.na(node_predator) == F) |> 
  left_join(biomasses, by = join_by(node_prey == node_name)) |> 
  group_by(node_predator) |> 
  mutate(rel_biomass = biomass / sum(biomass, na.rm = T)) |> 
  mutate(forage_ratio = ifelse(!is.na(a) & !is.na(h), (a * rel_biomass) / (1 + a * h * rel_biomass) / (rel_biomass), average_forage_ratio),
         forage_ratio = ifelse(is.na(forage_ratio) == T,
                               0, forage_ratio),
         # Check why NAs are introduced here!
         weight = (rel_biomass * forage_ratio) / sum(rel_biomass * forage_ratio, na.rm = T)) |> 
  ungroup() |> 
  
  # Make Table Graph
  select(node_predator, node_prey, forage_ratio, weight) |> 
  as_tbl_graph() |>
  activate(edges) |> 
  mutate(weight = replace_na(weight, 0)) |> 
  
  # Add node data  
  activate(nodes) |> 
  left_join(biomasses, by = join_by(name == node_name)) |> 
  mutate(biomass = replace_na(biomass, 0)) |> 
  left_join(bodymasses, by = join_by(name == node_name)) |> 
  mutate(bodymass = replace_na(bodymass, 1)) |> 
  left_join(node_data, by = join_by(name == node_name)) |> 
  
  # Calculate losses
  # TODO <- This calculation may not be correct!
  mutate(losses = exp(slope * log(bodymass) + intercept - tkonst),
         losses = ifelse(is.infinite(losses), 0, losses))

graph
# Tidygraph to fluxweb ----
# This is the end goal - not there yet.

fluxes <- 
  fluxing(mat = as_adjacency_matrix(graph, attr = "weight", sparse = FALSE),
          biomasses = pull(graph, biomass),
          losses = pull(graph, losses),
          efficiencies = pull(graph, efficiencies),
          bioms.prefs = FALSE,
          ef.level = "pred",
          bioms.losses = TRUE) |> 
  # Add biomass to node data again
  as_tbl_graph() |> activate(nodes) |> 
  activate(nodes) |> 
  left_join(biomasses, by = join_by(name == node_name)) |> 
  mutate(biomass = replace_na(biomass, 0))

# Principal visualization:
fluxes |> 
  ggraph() +
  geom_edge_link(aes(width = weight), arrow = arrow(length = unit(3, 'mm')),
                 alpha = 0.2) +
  geom_node_point(aes(size = biomass)) +
  geom_node_text(aes(label = name), angle = -90, size = 2.5, nudge_y = -0.2) +
  theme_graph() +
  #scale_color_manual(values = colors) +
  theme(legend.position = "none")




