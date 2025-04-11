#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))


#Prepare data ----
# Check if biomass data exists, otherwise create it.
if(!file.exists(file.path("data", "processed", "interpolation", "weekly_biomasses.csv"))) {
  system(paste("nohup Rscript", file.path("code", "InterpolateWeekly.R")))
}

# Check if predator selectivity data exists, otherwise create it.
if(!file.exists(file.path("data", "processed", "predator_selectivity.csv"))) {
  system(paste("nohup Rscript", file.path("code", "CombineMetabarcoding.R")))
}


# Merge selectivity and biomass ----
biomass <-
  read_csv(file.path("data", "processed", "interpolation", "weekly_biomasses.csv")) |> 
  group_by(sample_week, station_name) |>
  mutate(rel_biomass = biomass / sum(biomass)) |> 
  select(node_prey = node_name, sample_week, station_name, biomass)
ForageRatios <- 
  read_csv(file.path("data", "processed", "predator_selectivity.csv")) |> 
  ungroup() |>
  mutate(rra_gut = replace_na(rra_gut, 0),
         rra_env = replace_na(rra_env, 0),
         ForageRatio = rra_gut/rra_env) |>
  filter(rra_env > 0) |> 
  left_join(biomass) |> 
  group_by(sample_id) |> 
  mutate(rel_biomass = biomass / sum(biomass))

# Model ForageResponse ----