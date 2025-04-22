#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))

# The goal is to have a dataset that looks like:
# sample_id | node_predator | node_prey | barcode | predator_type | sample_date | sample_week | location | station_name | rra_gut | rra_env

# Fish forage ratio ------------------------------------------------------------
## COI -------------------------------------------------------------------------
df_COI <-
  # COI metadata
  read_csv(file.path("data", "raw", "fish_coi_metadata.csv"), show_col_types = FALSE) |>
  # The column sample_name contains multiple information
  extract(sample_name, into = c("station_ID", "sample_ID", "barcode"), 
          regex = "([^.]*)\\.(.*)_(.*)", remove = FALSE) |>
  full_join(# Combine ASV table with metadata
    read_csv(file.path("data","raw", "fish_coi_asv.csv"), show_col_types = FALSE) |>
      pivot_longer(cols = where(is.numeric),  # Pivot  numeric columns (abundance values)
                   names_to = "source_material_ID", 
                   values_to = "Abundance"),
    by = "source_material_ID"
  ) |>
  # Remove host reads
  filter(
    !(organism == "Sprattus sprattus" & Genus == "Sprattus") &
      !(organism == "Clupea harengus" & Genus == "Clupea") &
      !(organism == "Gasterosteus aculeatus" & Family == "Gasterosteidae")
  ) |>
  #  Keep genus that are in the node_name dataset
  mutate(Taxa = case_when(
    Genus %in% c("Acartia", "Eurytemora", "Evadne", "Pseudocalanus", "Synchaeta", "Temora") ~ Genus,
    TRUE ~ "Other"  # Default case
  )) |> filter(Taxa != "Other") |> 
  # remove sample with less than 10 000 reads
  group_by(sample_ID) |>
  filter(sum(Abundance, na.rm = T) > 10000) |> 
  # work in relative abundance
  mutate(RRA = Abundance / sum(Abundance, na.rm = T)) |>
  ungroup()

## 18S -------------------------------------------------------------------------
df_18S <-
  read_csv(file.path("data", "raw", "fish_18s_metadata.csv"), show_col_types = FALSE) |>
  # The column sample_name contains multiple information
  extract(sample_name, into = c("station_ID", "sample_ID", "barcode"), 
          regex = "([^.]*)\\.(.*)_(.*)", remove = FALSE) |> 
  full_join(# Combine ASV table with metadata
    read_csv(file.path("data", "raw", "fish_18s_asv.csv"), show_col_types = FALSE) |>
      pivot_longer(cols = where(is.numeric),  # Pivot  numeric columns (abundance values)
                   names_to = "source_material_ID", 
                   values_to = "Abundance"),
    by = "source_material_ID"
  ) |>
  #remove host reads
  filter(!(organism != "WP2" & Family == "Teleostei")) |> 
  # Keep genus that are in the node_name dataset
  mutate(Taxa = case_when(
    Genus %in% c("Acartia", "Centropages", "Eurytemora", "Evadne", "Pseudocalanus", "Synchaeta", "Temora") ~ Genus,
    Family == "Rotifera_XX" ~ "Synchaeta",  # Most Rotifera_XX are Synchaeta except specific ASVs
    ASV_ID %in% c("94ebe204e11f85a07820a7efef4f32d5", "d0cd100a81e83b5d909b9001096f202f") ~ "Other",
    TRUE ~ "Other"  # Default case
  )) |> filter(Taxa != "Other") |> 
  # remove sample with less than 10000 reads
  group_by(sample_ID) |>
  filter(sum(Abundance, na.rm = T) > 10000) |> 
  # work in relative abundance
  mutate(RRA = Abundance / sum(Abundance, na.rm = T)) |>
  ungroup() 
## 16S -------------------------------------------------------------------------
df_16S <-
  read_csv(file.path("data", "raw", "zooplankton_16s_metadata.csv"), show_col_types = FALSE) |>
  full_join(# Combine ASV table with metadata
    read_csv(file.path("data", "raw", "zooplankton_16s_asv.csv"), show_col_types = FALSE) |>
      pivot_longer(cols = where(is.numeric),  # Pivot  numeric columns (abundance values)
                   names_to = "library_ID", 
                   values_to = "Abundance"),
    by = "library_ID"
  ) |> 
  group_by(library_ID) |>
  filter(sum(Abundance, na.rm = T) > 10000) |> 
  ungroup() |> 
  # Keep genus that are in the node_name dataset
  mutate(Taxa = case_when(
    Order %in% c("Chaetocerotales", "Peridiniales", "Pyrenomonadales", "Thalassiosirales", "Eutreptiales", "Prymnesiales", "Pyramimonadales", "Cryptomonadales", "Bacillariales") ~ Order,
    Species == "Nodularia_PCC-9350" ~ "Nodularia",  # Most Rotifera_XX are Synchaeta except specific ASVs
    Species %in% c("Nostocaceaex", "Aphanizomenon_NIES81", "Aphanizomenon_MDT14a") ~ "Aphanizomenonaceae",
    Species %in% c("Cyanobium_PCC-6307", "Synechococcus_CC9902") ~ "Cyanobiaceae",
    Species == "Pseudanabaena_PCC-7429" ~ "Pseudanabaena",
    TRUE ~ "Other"  # Default case
  )) |>
  filter(Taxa != "Other",
         organism %in% c("Acartia", "Bosmina", "Centropages", "Eurytemora", "Evadne", "Pseudocalanus", "Seawater", "Synchaeta_baltica", "Temora")) |> 
  group_by(library_ID) |>
  # work in relative abundance
  mutate(RRA = Abundance / sum(Abundance, na.rm = T)) |>
  ungroup() 
## Combined --------------------------------------------------------------------
fish_zooplankton <-
  df_COI |>
  mutate(station_name = "BY31 LANDSORTSDJ",
         organism = case_when(
           organism == "Clupea harengus" ~ "Clupea",
           organism == "Sprattus sprattus" ~ "Sprattus",
           organism == "Gasterosteus aculeatus" ~ "Gasterosteus",
           organism == "WP2" ~ "environment"
         ),
         depth = NA) |>
  select(sample_ID, collection_date, station_ID, station_name, depth, barcode, sequence, organism, Taxa, RRA) |> 
  rename("sample_id" = sample_ID,
         "trawl_id" = station_ID,
         "ASV" = sequence,
         "node_predator" = organism,
         "node_prey" = Taxa,
         "rra" = RRA) |> 
  mutate(rra = replace_na(rra, 0)) |> 
  rbind(
    df_18S |> 
      mutate(station_name = "BY31 LANDSORTSDJ",
             organism = case_when(
               organism == "Clupea harengus" ~ "Clupea",
               organism == "Sprattus sprattus" ~ "Sprattus",
               organism == "Gasterosteus aculeatus" ~ "Gasterosteus",
               organism == "WP2" ~ "environment"
               ),
             depth = NA) |>
      select(sample_ID, collection_date, station_ID, station_name,depth, barcode, sequence, organism, Taxa, RRA) |> 
      rename("sample_id" = sample_ID,
             "trawl_id" = station_ID,
             "ASV" = sequence,
             "node_predator" = organism,
             "node_prey" = Taxa,
             "rra" = RRA) |> 
      mutate(rra = replace_na(rra, 0)))
zooplankton_phytoplankton <-
  df_16S |>
  filter(geo_loc_name %in% c("By31", "By2", "By5", "By16", "By15")) |> 
  mutate(
    station_name = case_when(
      geo_loc_name == "By31" ~ "BY31 LANDSORTSDJ",
      geo_loc_name == "By2" ~ "BY2 ARKONA",
      geo_loc_name == "By5" ~ "BY5 BORNHOLMSDJ",
      geo_loc_name %in% c("By16", "By15") ~ "BY15 GOTLANDSDJ"
      ),
    organism = case_when(
      organism == "Synchaeta_baltica" ~ "Synchaeta",
      organism == "Seawater" ~ "environment",
      organism %in% c("Eurytemora", "Evadne", "Temora", "Bosmina", "Pseudocalanus", "Centropages","Acartia") ~ organism
    ),
    barcode = "16S") |> 
  select(library_ID, collection_date,geo_loc_name, station_name, depth, barcode, ASV, organism, Taxa, RRA) |>
  rename("sample_id" = library_ID,
         "trawl_id" = geo_loc_name,
         "node_predator" = organism,
         "node_prey" = Taxa,
         "rra" = RRA) |> 
  mutate(rra = replace_na(rra, 0))

full_combined <-
  fish_zooplankton |> 
  rbind(zooplankton_phytoplankton) |> 
  mutate(sample_week = floor_date(collection_date, unit = "week", week_start = 1))
# split into predator and environment datasets
predator_df <-
  full_combined |> 
  filter(node_predator != "environment")
environment_df <-
  full_combined |> 
  filter(node_predator == "environment") |>
  group_by(collection_date,sample_week, trawl_id, station_name, barcode, ASV, node_prey, depth) |> 
  summarise(rra_tmp = mean(rra, na.rm = T), .groups = "drop_last") |> 
  summarise(rra_env = mean(rra_tmp, na.rm = T), .groups = "drop")
# Combine both datasets
predator_selectivity <-
  predator_df |> 
  left_join(environment_df,
            by = c("collection_date","sample_week", "trawl_id", "station_name", "barcode", "ASV", "node_prey"),
            relationship = "many-to-many") |>
  rename("rra_gut" = rra) |> 
  select(sample_id, collection_date, sample_week, station_name, barcode, ASV, node_predator, node_prey, rra_gut, rra_env) |> 
  distinct() |>
  filter(rra_env > 0) |> 
  group_by(sample_id,barcode) |> 
  mutate(rra_gut = rra_gut / sum(rra_gut, na.rm = T),
         rra_env = rra_env / sum(rra_env, na.rm = T)) |> 
  group_by(sample_id, collection_date, sample_week, station_name,node_predator,node_prey, barcode) |> 
  summarise(rra_gut = sum(rra_gut, na.rm = T),
            rra_env = sum(rra_env, na.rm = T),
            .groups = "drop_last") |> 
  summarise(rra_gut = mean(rra_gut, na.rm = T),
            rra_env = mean(rra_env, na.rm = T),
            .groups = "drop")

predator_selectivity |> 
  write_csv(file = file.path("data", "processed", "predator_selectivity.csv"))
