#!/usr/bin/env Rscript

# Dowload the shapefile data if not downloaded yet -----------------------------
if (!dir.exists(file.path("data", "imported", "ICES_areas")) |
    !dir.exists(file.path("data", "imported", "ICES_rectangles"))) {
  
  # 1. ICES Areas
  # The URL for getting the data for ICES areas
  areas_url <- "https://gis.ices.dk/shapefiles/ICES_areas.zip"
  # and the path to save it
  areas_zip_path <- file.path("data", "imported", "ICES_areas.zip")
  
  # Download
  download.file(areas_url, areas_zip_path, mode = "wb")
  
  # Unzip and delete the zip file
  areas_unziped_dir <- file.path("data", "imported", "ICES_areas")
  unzip(areas_zip_path, exdir = areas_unziped_dir)
  unlink(areas_zip_path)
  
  rm(areas_url, areas_zip_path, areas_unziped_dir)
  
  # 2. ICES rectangles
  # The URL for getting the data for ICES areas
  rect_url <- "https://gis.ices.dk/shapefiles/ICES_rectangles.zip"
  # and the path to save it
  rect_zip_path <- file.path("data", "imported", "ICES_rectangles.zip")
  
  # Download
  download.file(rect_url, rect_zip_path, mode = "wb")
  
  # Unzip and delete the zip file
  rect_unziped_dir <- file.path("data", "imported", "ICES_rectangles")
  unzip(rect_zip_path, exdir = rect_unziped_dir)
  unlink(rect_zip_path)
  
  rm(rect_url, rect_zip_path, rect_unziped_dir)
}

# Prepare the data -------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(sf))

# Names of the ices rectangle
ices_rect <- read_csv(file.path("data", "raw", "fish_parameters.csv"), show_col_types = FALSE) |> 
  pull(ICES_rect) |>
  unique()

# Coordinates of the SMHI stations
smhi_stations <-
  read_csv(file.path("data", "processed", "shark", "phytoplankton.csv"), show_col_types = FALSE) |>
  select(station_name, sample_latitude_dd, sample_longitude_dd) |>
  unique() |>
  slice(1:4)

# Coordinates where fish gut were estimated
spras <-
  read_csv(file.path("data", "raw", "fish_coi_metadata.csv"), show_col_types = FALSE) |> 
  rbind(read_csv(file.path("data", "raw", "fish_18s_metadata.csv"), show_col_types = FALSE)) |>
  filter(month(collection_date) == 05,
         collection_method == "Trawl") |> 
  select(lat_lon) |>
  separate(lat_lon, into = c("sample_latitude_dd", "sample_longitude_dd"), sep = " ") |> 
  mutate(across(c(sample_latitude_dd, sample_longitude_dd), as.numeric),
         sample_latitude_dd = round(sample_latitude_dd, 4),
         sample_longitude_dd = round(sample_longitude_dd, 4)) |> 
  unique()

# Baltic Sea area
baltic_sea_shp <-
  read_sf(list.files(file.path("data", "imported", "ICES_areas"), pattern = "\\.shp$", full.names = TRUE)) |> 
  filter(SubDivisio %in% 24:32) |>
  group_by(Major_FA) |> 
  summarise(geometry = st_union(geometry)) |> 
  ungroup()
rectangle_shp <- # Shapefile file with ices rectangle
  read_sf(list.files(file.path("data", "imported", "ICES_rectangles"), pattern = "\\.shp$", full.names = TRUE)) |> 
  filter(ICESNAME %in% ices_rect)

# Combine the ices statistical rectangle to the Baltic Sea shapefile to get the fish sampling area
if (st_crs(baltic_sea_shp) != st_crs(rectangle_shp)) {
  rectangle_shp <- st_transform(rectangle_shp, st_crs(baltic_sea_shp))
}
fish_shp <- # Combined dataset
  st_intersection(baltic_sea_shp, rectangle_shp)

# Transform the SMHI stations coordinates in a shapefile format
smhi_shp <-
  smhi_stations |> 
  st_as_sf(coords = c("sample_longitude_dd", "sample_latitude_dd"),
           crs = 4326)

# Transform the coordinates where fish gut were estimated in a shapefile format
spras_shp <-
  spras |> 
  st_as_sf(coords = c("sample_longitude_dd", "sample_latitude_dd"),
           crs = 4326)

# Plot the Baltic Sea map ------------------------------------------------------
suppressPackageStartupMessages(library(ggsflabel))

map <-
  ggplot() +
  # Start with the Baltic Sea contours
  geom_sf(data = baltic_sea_shp, fill = "white", color = "black") +
  # Add where fish biomass was estimated
  geom_sf(data = fish_shp,
          mapping = aes(fill = "Fish biomass"),
          color = "black") +
  # Add plankton sampling
  ggsflabel::geom_sf_text_repel(data = smhi_shp,
                                 mapping = aes(label = station_name),
                                 seed = 10, force = 100,
                                col = "#235452") +
  geom_sf(data = smhi_shp,
          mapping = aes(fill = "Plankton sampling"),
          shape = 21, size = 3) +
  # Add fish sampling
  geom_sf(data = spras_shp,
          mapping = aes(fill = "Fish sampling"),
          shape = 21, size = 2) +
  # Manual fill scale with labels and colors
  scale_fill_manual(
    name = NULL,
    values = c("Fish biomass" = "#F49E4C",
               "Plankton sampling" = "#235452",
               "Fish sampling" = "#F2EFC9")) +
  # Custom legend symbols
  guides(fill = guide_legend(
    override.aes = list(
      shape = c(NA, 21, 21),
      size = c(4, 2, 3),
      color = "black"))) +
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_text(color = "black", size = 10),
        panel.grid.major = element_line(color = "black",
                                        linewidth = 0.2))
# Save the map -----------------------------------------------------------------
## Create necessary directories if they do not exist (output)
if (!dir.exists(file.path("output"))) {
  dir.create(file.path("output"))
}
if (!dir.exists(file.path("output", "figure"))) {
  dir.create(file.path("output", "figure"))
}

ggsave(plot = map,
       filename = file.path("output", "figure", "map.pdf"),
       width = 6.5, height = 7.25)