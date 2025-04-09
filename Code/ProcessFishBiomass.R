#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))


read_csv(file.path("Data", "raw_fish_parameters.csv")) |> 
  group_by(node_name, year) |> 
  summarise(biomass = mean(biomass, na.rm = T),
            bodymass = mean(bodymass, na.rm = T),
            .groups = "drop") |> 
  full_join(
    tibble(
      sample_week = read_csv(file.path("Import", "SharkWeb","weekly_zooplankton.csv")) |> pull(sample_week) |> unique()) |> 
      mutate(year = year(sample_week)),
    by = "year",
    relationship = "many-to-many") |> 
  write_csv(file = file.path("Data", "weekly_fish_parameters.csv"))

