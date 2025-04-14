#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))

# Select Station and Day ----

sel_sample_week <- date("2007-03-12")
sel_station_name <- "BY31 LANDSORTSDJ"

# Calculate Weights ----

# Subset prey biomass
biomasses <-
  read_csv(file = file.path("data", "processed", "interpolation", "weekly_biomasses.csv")) |> 
  filter(sample_week == sel_sample_week,
         station_name == sel_station_name) |> 
  select(node_prey = node_name, biomass_prey = biomass)

# Calculate forage ratio and weights
# (Without bootstrap variation here)
weights <-
  read_csv(file = file.path("data", "processed", "forage_ratio.csv")) |>
  filter(is.na(node_predator) == F) |> 
  left_join(biomasses) |> 
  group_by(node_predator) |> 
  mutate(rel_biomass_prey = biomass_prey / sum(biomass_prey, na.rm = T)) |>
  mutate(ForageRatio = (a * rel_biomass_prey) / (1 + a * h * rel_biomass_prey) / (rel_biomass_prey),
         ForageRatio = ifelse(is.na(ForageRatio) == T,
                              average_forage_ratio, ForageRatio),
         # Check why NAs are introduced here!
         W = (rel_biomass_prey * ForageRatio) / sum(rel_biomass_prey * ForageRatio, na.rm = T)) |>
  ungroup()

# Turn Into Adjacency matrix
weights |> 
  select(node_predator, node_prey, W) |> 
  pivot_wider(names_from = node_predator, values_from = W, values_fill = 0) |> 
  column_to_rownames("node_prey") |> 
  as.matrix()



weights |> 
  ggplot() +
  geom_bar(aes(node_predator, W, fill = node_prey), stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
