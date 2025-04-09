#!/usr/bin/env Rscript

cat("\nRunning InterpolateWeekly.R\n")

suppressPackageStartupMessages(library(tidyverse))

if (!dir.exists(file.path("data", "processed"))) {
  dir.create(file.path("data", "processed"))
}

if (!dir.exists(file.path("data", "processed", "interpolation"))) {
  dir.create(file.path("data", "processed", "interpolation"))
}

if (!dir.exists(file.path("data", "processed", "shark"))) {
  dir.create(file.path("data", "processed", "shark"))
}



#this will download the data if they have not been downloaded yet, or will update them if a new version was released (it takes 20-30 min): 
if(!file.exists(file.path("data", "processed", "shark", "phytoplankton.csv")) |
   !file.exists(file.path("data", "processed", "shark", "zooplankton.csv")) |
   !file.exists(file.path("data", "processed", "shark", "temperature.csv")) |
   !file.exists(file.path("data", "processed", "shark", "picoplankton.csv"))) {
  getmonitoring <- file.path("code", "GetMonitoringData.R")
  # Get the monitoring data
  system(paste("nohup Rscript", getmonitoring))
}

# Zooplankton ------------------------------------------------------------------
zooplankton <-
  # Load zooplankton data

  read_csv(file.path("data", "processed", "shark", "zooplankton.csv"), show_col_types = FALSE) |>
  mutate(sample_week = floor_date(sample_date, unit = "week", week_start = 1),
         Month = month(sample_week),
         # Assign seasons based on the month
         season = case_when(
           Month %in% 1:3 ~ "winter",
           Month %in% 4:6 ~ "spring", 
           Month %in% 7:9 ~ "summer",
           Month %in% 10:12 ~ "fall"
         )) |>   # Change date to the first day of the week
  # Apply the filters
  filter(
    # From 0 to 60 m depth
    sample_min_depth_m %in% c(0, 30),
    sample_max_depth_m %in% c(30, 60)
    ) |>
  
  mutate(
    value = value * (sample_max_depth_m - sample_min_depth_m),
    unit = "ind/m2" # <---- from ind/m3 to ind/m2
    ) |>
  # Join the bodymass dataset
  right_join(
    read_csv(file.path("data","raw", "zooplankton_bodymass.csv"), show_col_types = FALSE),
    by = c("station_name", "sex_code", "dev_stage_code", "taxon_genus", "taxon_species", "season")
    ) |> 
#> 
  # Calculate biomass (g/m²) using abundance and body mass
  mutate(biomass_g = value * bodymass) |>
  # Summarise the dataset so we have one value per week
  filter(!is.na(sample_week)) |>
  group_by(sample_week, taxon_genus, station_name, sample_date) |>
  # Filter out the samples that were not collected in two depth layers
  filter(n_distinct(sample_min_depth_m) == 2) |> 
  # Add together the values that were collected the same date but at two different depth strata
  summarise(abundance_value = sum(value, na.rm = T),
            biomass_value = sum(biomass_g, na.rm = T),
            .groups = "drop_last") |> 
  # Summarise to have only one value per week
  summarise(abundance_ind.m2 = mean(abundance_value, na.rm = T),
            biomass_g.m2 = mean(biomass_value, na.rm = T),
            .groups = "drop") |>
# Interpolating zooplankton data and performing necessary transformations
  # Reshape to wide format
  pivot_wider(names_from = taxon_genus, values_from = c(abundance_ind.m2, biomass_g.m2)) |> 
  
  # Generate a complete sequence of weekly sample dates and join
  complete(sample_week = seq.Date(min(sample_week), max(sample_week), by = "week")) |>
  
  # Arrange by date and interpolate missing values
  arrange(sample_week) |> 
  mutate(across(-c(sample_week, station_name), ~ zoo::na.approx(.x, na.rm = FALSE)),
         station_name = "BY31") |>
  
  # Reshape back to long format
  pivot_longer(cols = -c(sample_week, station_name), names_to = "parameter", values_to = "value") |>
  
  # Split "parameter" into "Parameter" and "Taxa"
  separate(parameter, into = c("parameter", "unit", "Taxa"), sep = "_") |>
  select(- unit) |> 
  # Add Year and Week columns by extracting them from 'sample_date'
  mutate(Year = year(sample_week),
         Week = isoweek(sample_week)) |>
  
  # Reshape back to wide format using 'Parameter' to create columns
  pivot_wider(names_from = parameter, values_from = value) |>
  
  # Calculate Bodymass per individual by dividing Biomass by Abundance
  mutate(bodymass = biomass / abundance) |> 
  na.omit() |> # <- remove rows with NA 
  rename("node_name" = Taxa,
         "year" = Year)
zooplankton |>
  select(node_name, year, station_name, bodymass, sample_week) |> 
  write_csv(file.path("data", "processed", "interpolation", "zooplankton_bodymass.csv"))
zooplankton |>
  select(node_name, year, station_name, biomass, sample_week) |> 
  write_csv(file.path("data", "processed", "interpolation", "zooplankton_biomass.csv"))
cat("\nZooplankton weekly interpolated\n")


# Picophytoplankton ------------------------------------------------------------------
library(zoo)
# Read and filter
picoplankton <-
  read_csv(file.path("data", "processed", "shark", "picoplankton.csv")) |> 
  filter(scientific_name == "Synechococcus") |> 
  group_by(sample_id, station_name, sample_date) |> 
  summarise(value = sum(value)) |>
  filter(sample_date < date("2021-01-01"),
         sample_date > date("2018-10-30")) |> 
  arrange(sample_date) |> 
  group_by(station_name) |> 
  complete(sample_date = seq.Date(min(sample_date), max(sample_date), by = "week")) |>  # Fill in missing weekly dates
  arrange(sample_date) |> 
  mutate(value = na.approx(value, sample_date, na.rm = FALSE)) |>   # Linear interpolation
  ungroup() |> 
  filter(is.na(value) == F) |> 
  mutate(week = isoweek(sample_date)) |>   # Extract week of the year
  group_by(week) %>%
  summarize(value = mean(value, na.rm = TRUE), .groups = "drop")  # Average across years

# Phytoplankton ------------------------------------------------------------------



# Temperature ------------------------------------------------------------------
temperature <-
  read_csv(file.path("data", "processed","shark", "temperature.csv"), show_col_types = FALSE) |>  # Load the dataset
  mutate(sample_week = floor_date(sample_date, unit = "week", week_start = 1)) |>   # Change date to the first day of the week
  # Some checks and filters
  filter(
    station_name == "BY31 LANDSORTSDJ",
    sample_min_depth_m == sample_max_depth_m, # First check that the the depth is fixed while taking the temperature
    sample_min_depth_m %in% seq(from = 0, to = 60, by = 10) # And only select depth strata from 0 to 60m with 10m interval
    ) |>
  # Arrange the data
  group_by(sample_week, station_name, sample_min_depth_m) |> 
  summarise(value = mean(value, na.rm = T), .groups = "drop_last") |> # If there is 2 temperature record the same week, take the average value
  filter(n_distinct(sample_min_depth_m) == length(seq(from = 0, to = 60, by = 10))) |> # Make sure that all samples included in the analyses have all the temperature from 0 to 60m depth
  summarise(temperature = mean(value, na.rm = T), .groups = "drop") |> # and then take the average value of temperature from 0 to 60m depth

  # interpolation
  complete(sample_week = seq.Date(min(sample_week), max(sample_week), by = "week")) |> 
  # Arrange the data by sample_date to ensure proper chronological order
  arrange(sample_week) |>
  # Apply linear interpolation to all columns except 'sample_date'
  mutate(temperature = zoo::na.approx(temperature, na.rm = FALSE), 
         station_name = "BY31")
temperature |>
  write_csv(file.path("data","processed", "interpolation","temperature.csv"))

cat("\nTemperature weekly interpolated\n")

# Fish ----
## Biomass ----
read_csv(file.path("data", "raw","fish_parameters.csv")) |> 
  group_by(node_name, year) |> 
  summarise(biomass = mean(biomass, na.rm = T),
            bodymass = mean(bodymass, na.rm = T),
            .groups = "drop") |> 
  full_join(
    tibble(
      sample_week = read_csv(file.path("data", "processed", "interpolation", "zooplankton.csv")) |> pull(sample_week) |> unique()) |> 
      mutate(year = year(sample_week)),
    by = "year",
    relationship = "many-to-many") |>
  select(node_name, year, biomass, sample_week) |>
  write_csv(file = file.path("data", "processed", "interpolation", "fish_biomass.csv"))
## Bodymass ----
read_csv(file.path("data", "raw","fish_parameters.csv")) |> 
  group_by(node_name, year) |> 
  summarise(biomass = mean(biomass, na.rm = T),
            bodymass = mean(bodymass, na.rm = T),
            .groups = "drop") |> 
  full_join(
    tibble(
      sample_week = read_csv(file.path("data", "processed", "interpolation", "zooplankton.csv")) |> pull(sample_week) |> unique()) |> 
      mutate(year = year(sample_week)),
    by = "year",
    relationship = "many-to-many") |>
  select(node_name, year, bodymass, sample_week) |> 
  write_csv(file = file.path("data", "processed", "interpolation", "fish_bodymass.csv"))
