#!/usr/bin/env Rscript

# Create the directory for analyses
if (!dir.exists(file.path("data", "analyses")))
  dir.create(file.path("data", "analyses"))

## Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(tidygraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(fluxweb))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(abind))
suppressPackageStartupMessages(library(progressr))
# Check if bootstrap_forage_ratio exist, otherwise generate it.
if (
  !file.exists(file.path("data", "processed", "bootstrap_forage_ratio.csv"))
) {
  system(paste("nohup Rscript", file.path("code", "ModelForageResponse.R")))
}
cat("\nRunning ProcessFluxes.R\n")
# Source all needed functions for the analyses
source("./code/CalculateFluxes.R")

# Import the data
temperature <- read_csv(
  file = file.path("data", "processed", "interpolation", "temperature.csv"),
  show_col_types = FALSE
)
weekly_biomasses <- read_csv(
  file = file.path(
    "data",
    "processed",
    "interpolation",
    "weekly_biomasses.csv"
  ),
  show_col_types = FALSE
)
weekly_bodymass <- read_csv(
  file = file.path("data", "processed", "interpolation", "weekly_bodymass.csv"),
  show_col_types = FALSE
)
node_data <- read_csv(
  file = file.path("data", "raw", "node_data.csv"),
  show_col_types = FALSE
)
forage_ratio <- read_csv(
  file = file.path("data", "processed", "forage_ratio.csv"),
  show_col_types = FALSE
)
bootstrap_forage_ratio <- read_csv(
  file = file.path("data", "processed", "bootstrap_forage_ratio.csv"),
  show_col_types = FALSE
)

# Define the stations and the dates we want to compute the fluxes
station = "BY31 LANDSORTSDJ"
dates <- weekly_biomasses |>
  # Filter weekly_biomasses to only include samples from 2007 to end of 2023
  filter(
    sample_week > as.Date("2007-03-05"),
    sample_week <= as.Date("2023-12-04")
  ) |>
  pull(sample_week) |>
  unique()
# Null model ----
plan(multisession)
presence_absence_model <- tibble(sample_week = dates) |>
  mutate(
    flux_df = future_map(
      sample_week,
      ~ {
        flux_mat <- possibly(dna2flux, otherwise = NULL)(
          forage_ratio = forage_ratio,
          node_data = node_data,
          weekly_biomasses = weekly_biomasses,
          weekly_bodymass = weekly_bodymass,
          temperature = temperature,
          date = .x,
          station = station,
          as_graph = FALSE,
          presence_absence = TRUE,
          population_growth = FALSE
        )
        if (!is.null(flux_mat)) {
          flux_mat |>
            as.data.frame() |>
            rownames_to_column(var = "predator") |>
            pivot_longer(-predator, names_to = "prey", values_to = "flux") |>
            filter(flux > 0)
        }
      },
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  unnest(flux_df) |>
  mutate(station = station)

write_csv(
  presence_absence_model,
  file.path("data", "analyses", "timeseries_null.csv")
)


# Timeseries model ----
#
# Calculate fluxes using the forage ratio based on the 1000 iterations for each date...
# But it takes some time (about 2h per station), so it is better to do the computations once and save for later if we need to reuse

# Define the directory we want to save the bootrapped fluxes to be save
cache.dir = file.path("data", "analyses", "fluxes_array")
# Parallelize for faster computations
plan(multisession)

future_walk(
  dates,
  function(date) {
    cacheMyFluxes(
      cache.dir = cache.dir,
      bootstrap_forage_ratio = bootstrap_forage_ratio,
      node_data = node_data,
      weekly_biomasses = weekly_biomasses,
      weekly_bodymass = weekly_bodymass,
      temperature = temperature,
      date = date,
      station = station,
      as_graph = FALSE,
      presence_absence = FALSE,
      population_growth = FALSE
    )
  },
  .options = furrr_options(seed = TRUE),
  .progress = TRUE
)
cat(paste("\n All fluxes arrays are saved in", cache.dir))
# Create a tibble containing all dates and graph object with confidence interval and save it for later
# NOT SURE THAT IT IS USEFUL BUT LET'S SAVE IT ANYWAY FOR NOW
plan(multisession)
timeseries_fluxes <- tibble(
  sample_week = floor_date(date(dates), unit = "week", week_start = 1)
) |>
  # For each unique sampling week, compute the confidence graph in parallel
  mutate(
    conf_graph = future_map(
      sample_week,
      function(sample_date) {
        fluxingWithConfidence(
          bootstrap_forage_ratio = bootstrap_forage_ratio,
          node_data = node_data,
          weekly_biomasses = weekly_biomasses,
          weekly_bodymass = weekly_bodymass,
          temperature = temperature,
          date = sample_date,
          station = station,
          cache.dir = cache.dir,
          presence_absence = FALSE,
          population_growth = FALSE
        ) |>
          # Switch to node view (for adding metadata)
          activate(nodes) |>
          # Add sample week
          mutate(sample_week = sample_date)
      },
      .options = furrr_options(seed = TRUE) # Ensures reproducibility
    )
  )
# Save the daily_fluxes tibble as a rds file
write_rds(
  timeseries_fluxes,
  file = file.path("data", "analyses", "as_tbl_graph_timeseries_fluxes.rds")
)

# Modify the tibble to save timeseries_fluxes as a dataframe with sample_week, year, iso_week, predator, prey, mean +- CI fluxes
plan(multisession)
edges_df <- timeseries_fluxes |>
  pull(conf_graph) |>
  future_map_dfr(
    ~ .x |>
      activate(edges) |>
      mutate(across(c(flux_mean, flux_lower, flux_upper), ~ na_if(.x, 0))) |>
      extract_flux_long(),
    .options = furrr_options(seed = TRUE)
  ) |>
  mutate(
    iso_week = isoweek(sample_week),
    year = year(sample_week),
    station = station
  ) |>
  select(
    sample_week,
    year,
    iso_week,
    station,
    predator,
    prey,
    mean = flux_mean,
    upper = flux_upper,
    lower = flux_lower
  )
# Save as csv
write_csv(
  edges_df,
  file = file.path("data", "analyses", "timeseries_fluxes.csv")
)
# And get the node data too
node_df <- timeseries_fluxes |>
  pull(conf_graph) |>
  map_dfr(
    ~ .x |>
      activate(nodes) |>
      as_tibble()
  ) |>
  mutate(station_name = station)
write_csv(node_df, file = file.path("data", "analyses", "timeseries_nodes.csv"))
cat(paste("\n Timeseries fluxes for", station, "completed"))
# Time aggregated fluxes ----
# Create a tibble that have all the paths for each files
files <- tibble(path = list.files(cache.dir, full.names = TRUE)) |>
  mutate(
    filename = basename(path),
    date = str_extract(filename, "\\d{4}-\\d{2}-\\d{2}") |> as.Date(),
    year = year(date),
    iso_week = isoweek(date),
    station = str_extract(filename, "(?<=flux_)[^_]+")
  )
## Aggregate per iso_week ----
# Only select years that are after 2007
plan(multisession)
isoweek_fluxes <- files |>
  filter(year > 2007) |>
  group_nest(iso_week, station) |>
  mutate(
    flux = future_map(
      data,
      ~ {
        aggregated_fluxes <- aggregateFluxes(.x$path)
        if (is.null(aggregated_fluxes)) return(NULL)
        aggregated_fluxes
      },
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  select(-data) |>
  unnest(flux) |>
  group_by(predator, prey) |>
  filter(sum(mean, na.rm = T) > 0 & sum(upper, na.rm = T) > 0) |>
  ungroup()
# Save as csv
write_csv(
  isoweek_fluxes,
  file = file.path("data", "analyses", "isoweek_fluxes.csv")
)
cat(paste("\n Aggregated fluxes for", station, "based on isoweek completed"))
## Aggregate per year ----
plan(multisession)
annual_fluxes <- files |>
  filter(year > 2007) |> # As 2007 starts with values in March and misses January and February
  group_nest(year, station) |>
  mutate(
    flux = future_map(
      data,
      ~ {
        aggregated_fluxes <- aggregateFluxes(.x$path)
        if (is.null(aggregated_fluxes)) return(NULL)
        aggregated_fluxes
      },
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  select(-data) |>
  unnest(flux) |>
  group_by(predator, prey) |>
  filter(sum(mean, na.rm = T) > 0 & sum(upper, na.rm = T) > 0) |>
  ungroup()
# Save as csv
write_csv(
  annual_fluxes,
  file = file.path("data", "analyses", "annual_fluxes.csv")
)
cat(paste("\n Yearly aggregated fluxes for", station, "completed"))
## Aggregate per station over the entire timeseries ----
plan(multisession)
station_fluxes <- files |>
  group_nest(station) |>
  mutate(
    flux = future_map(
      data,
      ~ {
        aggregated_fluxes <- aggregateFluxes(.x$path)
        if (is.null(aggregated_fluxes)) return(NULL)
        aggregated_fluxes
      },
      .options = furrr_options(seed = TRUE)
    )
  ) |>
  select(-data) |>
  unnest(flux) |>
  group_by(predator, prey) |>
  filter(sum(mean, na.rm = T) > 0 & sum(upper, na.rm = T) > 0) |>
  ungroup()
# Save as csv
write_csv(
  station_fluxes,
  file = file.path("data", "analyses", "station_fluxes.csv")
)
cat(paste("\n Entire timeseries aggregated fluxes for", station, "completed"))
