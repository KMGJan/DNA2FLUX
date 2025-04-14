#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tidygraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(fluxweb))
suppressPackageStartupMessages(library(ggraph))

# Select Station and Day ----

sel_sample_week <- date("2007-03-12")
sel_station_name <- "BY31 LANDSORTSDJ"

# Calculate Weights ----
# Subset biomasses for day and station
biomasses <-
  read_csv(file = file.path("data", "processed", "interpolation", "weekly_biomasses.csv")) |> 
  filter(sample_week == sel_sample_week,
         station_name == sel_station_name) |> 
  select(node_name, biomass)

# Calculate forage ratio and weights
# (Without bootstrap variation here)
graph <-
  read_csv(file = file.path("data", "processed", "forage_ratio.csv")) |>
  filter(is.na(node_predator) == F) |> 
  left_join(biomasses, by = join_by(node_prey == node_name)) |> 
  group_by(node_predator) |> 
  mutate(rel_biomass = biomass / sum(biomass, na.rm = T)) |> 
  mutate(forage_ratio = (a * rel_biomass) / (1 + a * h * rel_biomass) / (rel_biomass),
         forage_ratio = ifelse(is.na(forage_ratio) == T,
                               average_forage_ratio, forage_ratio),
         # Check why NAs are introduced here!
         weight = (rel_biomass * forage_ratio) / sum(rel_biomass * forage_ratio, na.rm = T)) |> 
  ungroup() |> 
  select(node_predator, node_prey, forage_ratio, weight) |> 
  as_tbl_graph() |>
  activate(edges) |> 
  mutate(weight = replace_na(weight, 0)) |> 
  
  # Add biomasses to node data
  activate(nodes) |> 
  left_join(biomasses, by = join_by(name == node_name)) |> 
  mutate(biomass = replace_na(biomass, 0))


# Modify Node Data ----

graph









# Tidygraph to fluxweb ----
# This is the end goal - not there yet.

fluxes <- 
  fluxing(biomasses = as_adjacency_matrix(graph, attr = "W", sparse = FALSE),
          losses = pull(graph, losses),
          efficiencies = pull(graph, efficiencies),
          bioms.prefs = FALSE,
          ef.level = "prey",
          bioms.losses = TRUE)


