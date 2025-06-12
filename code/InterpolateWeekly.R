#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(zoo))

# Check that the data are present before interpolation
#this will download the data if they have not been downloaded yet, or will update them if a new version was released (it takes 20-30 min):
if (
  !file.exists(file.path("data", "processed", "shark", "phytoplankton.csv")) |
    !file.exists(file.path("data", "processed", "shark", "zooplankton.csv")) |
    !file.exists(file.path("data", "processed", "shark", "temperature.csv")) |
    !file.exists(file.path("data", "processed", "shark", "picoplankton.csv"))
) {
  getmonitoring <- file.path("code", "GetMonitoringData.R")
  # Get the monitoring data
  system(paste("nohup Rscript", getmonitoring))
}

message("Running InterpolateWeekly.R")

## Add output directory data/processed/interpolation  ----

if (!dir.exists(file.path("data", "processed", "interpolation"))) {
  dir.create(file.path("data", "processed", "interpolation"))
}

# Zooplankton ------------------------------------------------------------------

zooplankton <-
  # Load zooplankton data

  read_csv(
    file.path("data", "processed", "shark", "zooplankton.csv"),
    show_col_types = FALSE
  ) |>
  mutate(
    sample_week = floor_date(sample_date, unit = "week", week_start = 1),
    Month = month(sample_week),
    # Assign seasons based on the month
    season = case_when(
      Month %in% 1:3 ~ "winter",
      Month %in% 4:6 ~ "spring",
      Month %in% 7:9 ~ "summer",
      Month %in% 10:12 ~ "fall"
    )
  ) |> # Change date to the first day of the week
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
    read_csv(
      file.path("data", "raw", "zooplankton_bodymass.csv"),
      show_col_types = FALSE
    ),
    by = c(
      "station_name",
      "sex_code",
      "dev_stage_code",
      "taxon_genus",
      "taxon_species",
      "season"
    )
  ) |>

  # Calculate biomass (g/m²) using abundance and body mass
  mutate(biomass_g = value * bodymass) |>
  # Summarise the dataset so we have one value per week
  filter(!is.na(sample_week)) |>
  group_by(sample_week, taxon_genus, station_name, sample_date) |>
  # Filter out the samples that were not collected in two depth layers
  filter(n_distinct(sample_min_depth_m) == 2) |>
  # Add together the values that were collected the same date but at two different depth strata
  summarise(
    abundance_value = sum(value, na.rm = T),
    biomass_value = sum(biomass_g, na.rm = T),
    .groups = "drop_last"
  ) |>
  # Summarise to have only one value per week
  summarise(
    abundance_ind.m2 = mean(abundance_value, na.rm = T),
    biomass_g.m2 = mean(biomass_value, na.rm = T),
    .groups = "drop"
  ) |>

  # Interpolating zooplankton data and performing necessary transformations
  # Reshape to wide format
  pivot_wider(
    names_from = taxon_genus,
    values_from = c(abundance_ind.m2, biomass_g.m2),
    values_fill = 0
  ) |>
  group_by(station_name) |>
  # Generate a complete sequence of weekly sample dates and join
  complete(
    sample_week = seq.Date(min(sample_week), max(sample_week), by = "week")
  ) |>

  # Arrange by date and interpolate missing values
  arrange(sample_week, .by_group = TRUE) |>
  mutate(across(-c(sample_week), ~ na.approx(.x, na.rm = FALSE))) |>
  ungroup() |>
  # Reshape back to long format
  pivot_longer(
    cols = -c(sample_week, station_name),
    names_to = "parameter",
    values_to = "value"
  ) |>

  # Split "parameter" into "Parameter" and "Taxa"
  separate(parameter, into = c("parameter", "unit", "Taxa"), sep = "_") |>
  select(-unit) |>
  # Add Year and Week columns by extracting them from 'sample_date'
  mutate(Year = year(sample_week), Week = isoweek(sample_week)) |>

  # Reshape back to wide format using 'Parameter' to create columns
  pivot_wider(names_from = parameter, values_from = value) |>

  # Calculate Bodymass per individual by dividing Biomass by Abundance, when biomass = 0, bodymass will be NA, so linearly interpolation is needed
  group_by(station_name, Taxa) |>
  arrange(sample_week, .by_group = TRUE) |>
  mutate(
    bodymass = biomass / abundance,
    bodymass = na.approx(bodymass, na.rm = FALSE)
  ) |>
  ungroup() |>
  # na.omit() |> # <- remove rows with NA
  rename("node_name" = Taxa, "year" = Year)

zooplankton |>
  select(node_name, year, station_name, bodymass, sample_week) |>
  write_csv(file.path(
    "data",
    "processed",
    "interpolation",
    "zooplankton_bodymass.csv"
  ))
zooplankton |>
  select(node_name, year, station_name, biomass, sample_week) |>
  write_csv(file.path(
    "data",
    "processed",
    "interpolation",
    "zooplankton_biomass.csv"
  ))
message("Zooplankton weekly interpolated")


# Picophytoplankton ------------------------------------------------------------------

# Picocyanobacteria are not included in the phytoplankton dataset, but has been
# counted at a handful of stations separately.

# Read and filter
picoplankton <-
  read_csv(
    file.path("data", "processed", "shark", "picoplankton.csv"),
    show_col_types = FALSE
  ) |>
  filter(scientific_name == "Synechococcus") |>
  group_by(sample_id, station_name, sample_date) |>
  summarise(
    # Convert ugC to biomass g/m2
    # * 10 for sampling depth;  * 4 for carbon to biomass; * 0.001 ug to g
    biomass = sum(value) * 10 * 4 * 0.001,
    .groups = "drop"
  ) |>
  filter(sample_date < date("2021-01-01"), sample_date > date("2018-10-30")) |>
  arrange(sample_date) |>
  group_by(station_name) |>
  complete(
    sample_date = seq.Date(min(sample_date), max(sample_date), by = "week")
  ) |> # Fill in missing weekly dates
  arrange(sample_date) |>
  mutate(biomass = na.approx(biomass, sample_date, na.rm = FALSE)) |> # Linear interpolation
  ungroup() |>
  filter(is.na(biomass) == F) |>
  mutate(week_number = isoweek(sample_date)) |> # Extract week of the year
  group_by(week_number) %>%
  summarize(Cyanobiaceae = mean(biomass, na.rm = TRUE), .groups = "drop") # Average across years


# Phytoplankton ------------------------------------------------------------------
## Read shark and modify taxa to node names ----
node_names <-
  read_csv(file.path("data", "raw", "node_data.csv"), show_col_types = FALSE) |>
  pivot_longer(3:5, names_to = "data_type", values_to = "name") |>
  filter(is.na(name) == F) |>
  separate_rows(name, sep = ";")

order_pp <- node_names |>
  filter(
    data_type == "dyntaxa_name",
    type == "phytoplankton",
    tax_level == "order"
  ) |>
  select(node_name, taxon_order = name) |>
  right_join(
    read_csv(
      file.path("data", "processed", "shark", "phytoplankton.csv"),
      show_col_types = FALSE
    ),
    by = "taxon_order"
  ) |>
  filter(is.na(node_name) == F)

genus_pp <- node_names |>
  filter(
    data_type == "dyntaxa_name",
    type == "phytoplankton",
    tax_level == "genus"
  ) |>
  select(node_name, taxon_genus = name) |>
  right_join(
    read_csv(
      file.path("data", "processed", "shark", "phytoplankton.csv"),
      show_col_types = FALSE
    ),
    by = "taxon_genus"
  ) |>
  filter(is.na(node_name) == F)
## Interpolation ----
bind_rows(genus_pp, order_pp) |>
  mutate(
    sample_week = floor_date(sample_date, unit = "week", week_start = 1),
    # Convert ugC to biomass g/m2
    # * sampling depth;  * 4 for carbon to biomass; * 0.001 ug to g
    value = value * (sample_max_depth_m - sample_min_depth_m) * 4 * 0.001
  ) |>
  group_by(
    station_name,
    sample_week,
    node_name,
    shark_sample_id_md5,
    sample_max_depth_m
  ) |>
  summarise(biomass = sum(value), .groups = "drop_last") |>
  summarise(biomass = mean(biomass), .groups = "drop_last") |>
  # Average value for each sample_week, station
  summarise(biomass = mean(biomass, na.rm = T), .groups = "drop") |>

  # Reshape to wide format and merge with picoplankton by week
  pivot_wider(names_from = node_name, values_from = biomass, values_fill = 0) |>

  # Interpolating zooplankton data and performing necessary transformations
  group_by(station_name) |>

  # Generate a complete sequence of weekly sample dates and join
  complete(
    sample_week = seq.Date(min(sample_week), max(sample_week), by = "week")
  ) |>

  # Arrange by date and interpolate missing values
  arrange(sample_week, .by_group = TRUE) |>
  mutate(across(-c(sample_week), ~ na.approx(.x, na.rm = FALSE))) |>

  # Reshape to wide format and merge with picoplankton by week
  mutate(week_number = isoweek(sample_week)) |>
  ungroup() |>
  left_join(picoplankton, by = "week_number") |>
  select(!week_number) |>

  # Reshape back to long format
  pivot_longer(
    cols = -c(sample_week, station_name),
    names_to = "node_name",
    values_to = "biomass"
  ) |>

  select(node_name, sample_week, station_name, biomass) |>

  write_csv(file.path(
    "data",
    "processed",
    "interpolation",
    "phytoplankton_biomass.csv"
  ))


# Temperature ------------------------------------------------------------------
temperature <-
  read_csv(
    file.path("data", "processed", "shark", "temperature.csv"),
    show_col_types = FALSE
  ) |> # Load the dataset
  mutate(
    sample_week = floor_date(sample_date, unit = "week", week_start = 1)
  ) |> # Change date to the first day of the week
  # Some checks and filters
  filter(
    sample_min_depth_m == sample_max_depth_m, # First check that the the depth is fixed while taking the temperature
    sample_min_depth_m %in% seq(from = 0, to = 60, by = 10) # And only select depth strata from 0 to 60m with 10m interval
  ) |>
  # Arrange the data
  group_by(station_name, sample_week, sample_min_depth_m) |>
  summarise(value = mean(value, na.rm = T), .groups = "drop_last") |> # If there is 2 temperature record the same week, take the average value
  # Make sure that all samples included in the analyses have all the temperature from 0 to 60m depth for all stations and between 0 and 40 for BY2
  filter(
    (station_name != "BY2 ARKONA" &
      n_distinct(sample_min_depth_m) ==
        length(seq(from = 0, to = 60, by = 10))) |
      (station_name == "BY2 ARKONA" &
        n_distinct(sample_min_depth_m) ==
          length(seq(from = 0, to = 40, by = 10)))
  ) |>
  summarise(temperature = mean(value, na.rm = T), .groups = "drop_last") |> # and then take the average value of temperature from 0 to 60m depth

  # interpolation
  complete(
    sample_week = seq.Date(min(sample_week), max(sample_week), by = "week")
  ) |>
  # Arrange the data by sample_date to ensure proper chronological order
  arrange(sample_week, .by_group = TRUE) |>
  # Apply linear interpolation to all columns except 'sample_date'
  mutate(temperature = na.approx(temperature, na.rm = FALSE))

temperature |>
  write_csv(file.path("data", "processed", "interpolation", "temperature.csv"))

message("Temperature weekly interpolated")


# Fish ----

#' Read and Combine Fish Biomass Data Files
#'
#' Reads multiple fish biomass CSV files (semicolon-separated, with "." as decimal mark), extracts the year from each filename, and appends it as a column in the resulting data.
#'
#' @param list_of_files A character vector of file paths to `.csv` files containing fish biomass data. Each file must contain a 4-digit year in the filename (e.g., `fishdata_2020.csv`).
#'
#' @return A single data frame combining all the input files, with an additional column `year` indicating the year extracted from each file name.
#'
#' @examples
#' \dontrun{
#' files <- list.files("data/raw/BIAS Stickleback", full.names = TRUE)
#' biomass_data <- readFishBiomass(files)
#' }
#'
readFishBiomass <- function(list_of_files) {
  data_list <- lapply(list_of_files, function(file) {
    df <- read.csv2(file, dec = ".")

    # Extract the year from the file name
    year <- sub(".*?(\\d{4})\\.csv$", "\\1", file)

    # Add the year as a new column
    df$year <- as.integer(year)

    return(df)
  }) |>
    bind_rows()
}

stickleback <-
  list.files(file.path("data", "raw", "BIAS Stickleback"), full.names = TRUE) |>
  readFishBiomass() |>
  rename(
    "ICES_rect" = ICES,
    "N" = TotNo,
    "biomass" = Biom,
    "bodymass" = meanW
  ) |>
  select(-SD) |>
  mutate(node_name = "Gasterosteus")

herring <-
  list.files(
    file.path("data", "raw", "BIAS Herr Sprat"),
    pattern = "^biohe",
    full.names = TRUE
  ) |>
  readFishBiomass() |>
  # So we have everything in grams
  mutate(
    N = N * 1e6, # Abundance was in million -> now in individual
    biom = biom * 1e6
  ) |>
  rename("ICES_rect" = RECT, "biomass" = biom, "bodymass" = meanW) |>
  mutate(node_name = "Clupea")

sprat <-
  list.files(
    file.path("data", "raw", "BIAS Herr Sprat"),
    pattern = "^biosp",
    full.names = TRUE
  ) |>
  readFishBiomass() |>
  # So we have everything in grams
  mutate(
    N = N * 1e6, # Abundance was in million -> now in individual
    biom = biom * 1e6
  ) |>
  rename("ICES_rect" = RECT, "biomass" = biom, "bodymass" = meanW) |>
  mutate(node_name = "Sprattus")

# Area of ICES rectangle
area_df <-
  data.frame(
    ICES_rect = c(
      "42G6",
      "42G7",
      "43G6",
      "43G7",
      "43G8",
      "44G6",
      "44G7",
      "44G8",
      "45G6",
      "45G7",
      "45G8",
      "46G6",
      "46G7",
      "46G8",
      "47G8",
      "48G8"
    ),
    Area_nm2 = c(
      266.0,
      986.9,
      269.8,
      913.8,
      106.1,
      200.9,
      960.5,
      456.6,
      72.9,
      908.7,
      947.2,
      38.9,
      452.6,
      884.8,
      264.3,
      53.8
    )
  ) |>
  mutate(Area_m2 = Area_nm2 * 1852^2) # 1 nautical mile = 1852 meters.

# Combine all three datasets and calculate the biomass in g/m2 and the abundance in ind/m2
fish <-
  stickleback |>
  bind_rows(herring) |>
  bind_rows(sprat) |>
  left_join(area_df, by = "ICES_rect") |>
  mutate(biomass = biomass / Area_m2, abundance = N / Area_m2) |>
  select(node_name, ICES_rect, abundance, biomass, bodymass, year) |>
  arrange(node_name, year, ICES_rect) |>
  mutate(
    station_name = ifelse(
      ICES_rect %in% c("43G9", "43H0"),
      "BY15 GOTLANDSDJ",
      ifelse(
        ICES_rect %in% c("39G5", "39G6"),
        "BY5 BORNHOLMSDJ",
        ifelse(ICES_rect %in% c("45G8", "46G8"), "BY31 LANDSORTSDJ", ICES_rect)
      )
    )
  )

# Interpolation
fish |>
  group_by(node_name, year, station_name) |>
  summarise(
    biomass = mean(biomass, na.rm = T),
    bodymass = mean(bodymass, na.rm = T),
    .groups = "drop"
  ) |>
  filter(
    station_name %in%
      c("BY31 LANDSORTSDJ", "BY5 BORNHOLMSDJ", "BY15 GOTLANDSDJ")
  ) |>
  full_join(
    tibble(
      sample_week = read_csv(
        file.path(
          "data",
          "processed",
          "interpolation",
          "zooplankton_biomass.csv"
        ),
        show_col_types = FALSE
      ) |>
        pull(sample_week) |>
        unique()
    ) |>
      mutate(year = year(sample_week)),
    by = "year",
    relationship = "many-to-many"
  ) |>
  select(node_name, year, biomass, sample_week, station_name) |>
  write_csv(
    file = file.path("data", "processed", "interpolation", "fish_biomass.csv")
  )

## Bodymass ----
fish |>
  group_by(node_name, year, station_name) |>
  summarise(
    biomass = mean(biomass, na.rm = T),
    bodymass = mean(bodymass, na.rm = T),
    .groups = "drop"
  ) |>
  filter(
    station_name %in%
      c("BY31 LANDSORTSDJ", "BY5 BORNHOLMSDJ", "BY15 GOTLANDSDJ")
  ) |>
  full_join(
    tibble(
      sample_week = read_csv(
        file.path(
          "data",
          "processed",
          "interpolation",
          "zooplankton_bodymass.csv"
        ),
        show_col_types = FALSE
      ) |>
        pull(sample_week) |>
        unique()
    ) |>
      mutate(year = year(sample_week)),
    by = "year",
    relationship = "many-to-many"
  ) |>
  select(node_name, year, bodymass, sample_week, station_name) |>
  write_csv(
    file = file.path("data", "processed", "interpolation", "fish_bodymass.csv")
  )


# Merge Biomasses ----------

read_csv(
  file.path("data", "processed", "interpolation", "fish_biomass.csv"),
  show_col_types = FALSE
) |>
  bind_rows(read_csv(
    file.path("data", "processed", "interpolation", "zooplankton_biomass.csv"),
    show_col_types = FALSE
  )) |>
  bind_rows(read_csv(
    file.path(
      "data",
      "processed",
      "interpolation",
      "phytoplankton_biomass.csv"
    ),
    show_col_types = FALSE
  )) |>
  write_csv(
    file = file.path(
      "data",
      "processed",
      "interpolation",
      "weekly_biomasses.csv"
    )
  )

# Merge Bodymass ----------

read_csv(
  file.path("data", "processed", "interpolation", "fish_bodymass.csv"),
  show_col_types = FALSE
) |>
  bind_rows(read_csv(
    file.path("data", "processed", "interpolation", "zooplankton_bodymass.csv"),
    show_col_types = FALSE
  )) |>
  write_csv(
    file = file.path(
      "data",
      "processed",
      "interpolation",
      "weekly_bodymass.csv"
    )
  )

# Clean the environment
rm(list = ls())
