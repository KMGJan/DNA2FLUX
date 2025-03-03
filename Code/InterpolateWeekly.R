#!/usr/bin/env Rscript

cat("\nRunning InterpolateWeekly.R\n")

suppressPackageStartupMessages(library(tidyverse))

# this will download the data if they have not been downloaded yet, or will update them if a new version was released (it takes 20-30 min): 
if(!file.exists(file.path("Import", "SharkWeb", "BY31", "phytoplankton.csv")) |
   !file.exists(file.path("Import", "SharkWeb", "BY31", "zooplankton.csv")) |
   !file.exists(file.path("Import", "SharkWeb", "BY31", "temperature.csv"))) {
  getmonitoring <- file.path("Code", "GetMonitoringData.R")
  # Get the monitoring data
  system(paste("nohup Rscript", getmonitoring))
}

# Zooplankton ------------------------------------------------------------------
zooplankton <-
  # Load zooplankton data
  read_csv(file.path("Import", "SharkWeb", "BY31", "zooplankton.csv"), show_col_types = FALSE) |>
  mutate(sample_week = floor_date(sample_date, unit = "week", week_start = 1)) |>   # Change date to the first day of the week
  # Apply the filters
  filter(
    # Keep only selected taxonomic classes
    taxon_class %in% c("Maxillopoda", "Branchiopoda", "Eurotatoria"),
    # That are not nauplii     
    dev_stage_code != "NP",
    # From 0 to 60 m depth
    sample_min_depth_m %in% c(0, 30),
    sample_max_depth_m %in% c(30, 60),
    taxon_genus %in% c(
      "Acartia", "Temora", "Centropages", "Eurytemora", "Pseudocalanus", # Copepoda
      "Keratella", "Synchaeta", # Rotifera
      "Evadne", "Podon", "Pleopis", "Bosmina") # Cladocera
    ) |>
  
  # Arrange the data taxonomy as the bodymass data
  mutate(
    # Standardize development stage categories
    dev_stage_code = ifelse(dev_stage_code %in% c("C3", "C4"), "C4", as.character(dev_stage_code)),
    # Standardize station names
    station_name = ifelse(station_name == "BY31 LANDSORTSDJ", "BY31", as.character(station_name)),
    # Extract month from sample date
    Month = month(sample_week),
    # Assign seasons based on the month
    season = case_when(
      Month %in% 1:3 ~ "winter",
      Month %in% 4:6 ~ "spring", 
      Month %in% 7:9 ~ "summer",
      Month %in% 10:12 ~ "fall"
      ),
    # Assign Taxa
    Taxa = ifelse(taxon_species %in% c("Acartia bifilosa", "Acartia longiremis", "Acartia tonsa",
                                       "Podon intermedius", "Podon leuckartii"),
                       as.character(taxon_species),
                       ifelse(taxon_genus == "Pleopis", "Podon", as.character(taxon_genus))),
    value = value * (sample_max_depth_m - sample_min_depth_m),
    unit = "ind/m2", # <---- from ind/m3 to ind/m2
    across(c(sex_code, dev_stage_code), # Apply the transformation to both 'sex_code' and 'dev_stage_code' columns
           ~ ifelse(. %in% c("M", "F", "AD", "C1", "C4", "JV"), ., "NS"))) |>  # Keep valid values, replace others with "NS"
  # Join the bodymass dataset
  right_join(
    read_csv(file.path("Data", "bodymass.csv"), show_col_types = FALSE),
    by = c("station_name", "sex_code", "dev_stage_code", "season", "Taxa")
    ) |> 
  # All the taxa should now be at the genus level
  mutate(
    Taxa = case_when(
      str_detect(Taxa, "^Acartia") ~ "Acartia",
      str_detect(Taxa, "^Podon") ~ "Podon",
      TRUE ~ as.character(Taxa)
      ), 
    # Calculate biomass (g/m²) using abundance and body mass
    biomass_g = value * bodymass * 1e-6) |> # Convert from micrograms to grams (g)
  # Summarise the dataset so we have one value per week
  filter(!is.na(sample_week)) |>
  group_by(sample_week, Taxa, station_name, sample_date) |>
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
  pivot_wider(names_from = Taxa, values_from = c(abundance_ind.m2, biomass_g.m2)) |> 
  
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
  
  # Rename columns for clarity
  rename("Abundance_ind.m2" = abundance,
         "Biomass_g.m2" = biomass) |>
  
  # Calculate Bodymass per individual by dividing Biomass by Abundance
  mutate(Bodymass_g.ind = Biomass_g.m2 / Abundance_ind.m2) |> 
  na.omit() # <- remove rows with NA
zooplankton |>
  write_csv(file.path("Import", "SharkWeb", "weekly_zooplankton.csv"))
cat("\nZooplankton weekly interpolated\n")
# Temperature ------------------------------------------------------------------
temperature <-
  read_csv(file.path("Import", "Sharkweb", "BY31", "temperature.csv"), show_col_types = FALSE) |>  # Load the dataset
  mutate(sample_week = floor_date(sample_date, unit = "week", week_start = 1)) |>   # Change date to the first day of the week
  # Some checks and filters
  filter(
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
  write_csv(file.path("Import", "SharkWeb", "weekly_temperature.csv"))

cat("\nTemperature weekly interpolated\n")