#!/usr/bin/env Rscript

cat("\nRunning GetMonitoringData.R\n")

## Create necessary directories if they do not exist
#if (!dir.exists(file.path("data","imported", "sharkweb"))) {
#  dir.create(file.path("data","imported", "sharkweb"))
#}
raw_dir <- file.path("data","imported","sharkweb")
#if (!dir.exists(raw_dir)) {
#  dir.create(raw_dir)
#}
#SHARKphysical <- file.path("code", "SHARKphysical.sh")
## Run the data download scripts in the background using `nohup`
#system(paste("nohup bash", SHARKphysical, raw_dir)) # Download physical-chemical data
#
#SHARKplankton <- file.path("code", "SHARKplankton.sh")
#system(paste("nohup bash", SHARKplankton, raw_dir)) # Download plankton data





cat(paste("\nAll data packages are downloaded in", raw_dir, "and up to date\n"))

# Load necessary libraries with suppressed messages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# Combine all CSV files into a single output directory
options(readr.show_progress = FALSE)
output_dir <- file.path("data", "processed", "shark")
cat(paste("\nAll data packages are being combined and will be stored in", output_dir, "when completed\n"))

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define column types for reading CSV files
col_type <- list(
  delivery_datatype = col_character(), check_status_sv = col_character(), data_checked_by_sv = col_character(), visit_year = col_double(), visit_month = col_double(), station_name = col_character(), reported_station_name = col_character(), sample_location_id = col_character(), station_id = col_character(), sample_project_name_sv = col_character(), sample_orderer_name_sv = col_character(), visit_id = col_character(), visit_date = col_date(format = ""), shark_sample_id_md5 = col_character(), sample_date = col_date(format = ""), sample_time = col_time(format = ""), sample_enddate = col_date(), sample_endtime = col_date(), sample_latitude_dm = col_character(),  sample_longitude_dm = col_character(), sample_latitude_dd = col_double(), sample_longitude_dd = col_double(), water_depth_m = col_double(), visit_comment = col_character(), sample_id = col_character(), sample_min_depth_m = col_double(), sample_max_depth_m = col_double(), sampling_laboratory_name_sv = col_character(), sample_comment = col_character(), scientific_name = col_character(), species_flag_code = col_character(), dyntaxa_id = col_character(), aphia_id = col_character(), parameter = col_character(), value = col_double(), unit = col_character(), quality_flag = col_character(), calc_by_dc = col_character(), sex_code = col_character(), dev_stage_code = col_character(), trophic_type_code = col_character(), size_class = col_character(), method_documentation = col_character(), method_reference_code = col_character(), variable_comment = col_character(), analytical_laboratory_name_sv = col_character(), location_water_category = col_character(), location_water_district = col_character(), location_svar_sea_area_name = col_character(), location_svar_sea_area_code = col_character(), location_type_area = col_character(), location_sea_basin = col_character(), location_helcom_ospar_area = col_character(), location_economic_zone = col_character(), location_county = col_character(), location_municipality = col_character(), station_viss_eu_id = col_character(), taxon_kingdom = col_character(), taxon_phylum = col_character(), taxon_class = col_character(), taxon_order = col_character(), taxon_family = col_character(), taxon_genus = col_character(), taxon_species = col_character(), taxon_hierarchy = col_character(), taxon_red_list_category = col_character(), reported_scientific_name = col_character(), reported_parameter = col_character(), reported_value = col_character(), reported_unit = col_character(), dataset_name = col_character(), dataset_file_name = col_character()
  )

# Process Phytoplankton Data
lapply(
  list.files(raw_dir, pattern = "^SHARK_Phytoplankton.*\\.csv$", full.names = TRUE),
  read_csv,
  show_col_types = FALSE,
  col_types = col_type,
  col_select = c(shark_sample_id_md5, delivery_datatype, visit_year, station_name, sample_project_name_sv, sample_orderer_name_sv,
                 sample_date, sample_latitude_dd, sample_longitude_dd, sample_min_depth_m, sample_max_depth_m,
                 sampling_laboratory_name_sv, scientific_name, dyntaxa_id, aphia_id, parameter, value, unit,
                 sex_code, dev_stage_code, trophic_type_code, taxon_kingdom, taxon_phylum, taxon_class, 
                 taxon_order, taxon_family, taxon_genus, taxon_species, scientific_name)
  ) |>
  bind_rows() |>
  filter(station_name %in% c( "BY31 LANDSORTSDJ",
                              "BY5 BORNHOLMSDJ",
                              "BY15 GOTLANDSDJ",
                              "BY2 ARKONA"), unit == "ugC/l") |>  # Filter for station and unit
  write_csv(file.path(output_dir, "phytoplankton.csv"))  # Save processed file


list.files(raw_dir, pattern = "^SHARK_Phytoplankton.*\\.csv$", full.names = TRUE)
cat("\nPhytoplankton dataset done\n")

# Process Zooplankton Data
lapply(
  list.files(raw_dir, pattern = "^SHARK_Zooplankton.*\\.csv$", full.names = TRUE),
  read_csv,
  show_col_types = FALSE,
  col_types = col_type,
  col_select = c(delivery_datatype, visit_year, station_name, sample_project_name_sv, sample_orderer_name_sv,
                 sample_date, sample_latitude_dd, sample_longitude_dd, sample_min_depth_m, sample_max_depth_m,
                 sampling_laboratory_name_sv, scientific_name, dyntaxa_id, aphia_id, parameter, value, unit,
                 sex_code, dev_stage_code, trophic_type_code, taxon_kingdom, taxon_phylum, taxon_class,
                 taxon_order, taxon_family, taxon_genus, taxon_species)
) |> 
  bind_rows() |> 
  filter(station_name %in% c( "BY31 LANDSORTSDJ",
                              "BY5 BORNHOLMSDJ",
                              "BY15 GOTLANDSDJ",
                              "BY2 ARKONA"), unit == "ind/m3") |>  # Filter for station and unit
  write_csv(file.path(output_dir, "zooplankton.csv"))  # Save processed file

cat("\nZooplankton dataset done\n")

# Process Temperature Data
lapply(
  list.files(raw_dir, pattern = "^SHARK_PhysicalChemical.*\\.csv$", full.names = TRUE),
  read_csv,
  show_col_types = FALSE,
  col_types = col_type,
  col_select = c(delivery_datatype, visit_year, station_name, sample_project_name_sv, sample_orderer_name_sv,
                 sample_date, sample_latitude_dd, sample_longitude_dd, sample_min_depth_m, sample_max_depth_m,
                 sampling_laboratory_name_sv, parameter, value, unit)
) |>
  bind_rows() |>
  filter(station_name %in% c( "BY31 LANDSORTSDJ",
                              "BY5 BORNHOLMSDJ",
                              "BY15 GOTLANDSDJ",
                              "BY2 ARKONA"), parameter == "Temperature CTD") |>  # Filter for temperature data
  write_csv(file.path(output_dir, "temperature.csv"))  # Save processed file

cat("\nTemperature dataset done\n")

# Process Picoplankton Data
lapply(
  list.files(raw_dir, pattern = "^SHARK_Picoplankton.*\\.csv$", full.names = TRUE),
  read_csv,
  show_col_types = FALSE,
  col_types = col_type,
  col_select = c(sample_id, delivery_datatype, visit_year, station_name, sample_project_name_sv, sample_orderer_name_sv,
                 sample_date, sample_latitude_dd, sample_longitude_dd, sample_min_depth_m, sample_max_depth_m,
                 sampling_laboratory_name_sv, scientific_name, dyntaxa_id, aphia_id, parameter, value, unit,
                 sex_code, dev_stage_code, trophic_type_code, taxon_kingdom, taxon_phylum, taxon_class, 
                 taxon_order, taxon_family, taxon_genus, taxon_species)
) |>
  bind_rows() |>
  filter(station_name %in% c( "BY31 LANDSORTSDJ",
                              "BY5 BORNHOLMSDJ",
                              "BY15 GOTLANDSDJ",
                              "BY2 ARKONA"), unit == "ugC/l") |>  # Filter for station and unit
  write_csv(file.path(output_dir, "picoplankton.csv"))  # Save processed file

cat("\nPicoplankton dataset done\n")
