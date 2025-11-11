#!/usr/bin/env Rscript

# Fri Aug 15 08:57:18 2025 ------------------------------

# First ensure that all needed packages are installed, or install them
source(file.path("code", "InstallPackages.R"))

# Ensure that all datasets were generated
processedFluxes <- c(
  "annual_fluxes.csv",
  "as_tbl_graph_timeseries_fluxes.rds",
  "isoweek_fluxes.csv",
  "station_fluxes.csv",
  "timeseries_fluxes.csv",
  "timeseries_nodes.csv",
  "timeseries_null.csv"
)
fluxPaths <- file.path("data", "analyses", processedFluxes)
if (!all(file.exists(fluxPaths))) {
  message("The needed dataset to run the analyses are not processed yet")
  system(paste("Rscript", file.path("code", "ProcessFluxes.R")))
}
rm(processedFluxes, fluxPaths)
# Create the directories output/figure and output/table
if (!dir.exists(file.path("output", "figure")))
  dir.create(file.path("output", "figure"), recursive = T)
if (!dir.exists(file.path("output", "residuals")))
  dir.create(file.path("output", "residuals"), recursive = T)
if (!dir.exists(file.path("output", "table")))
  dir.create(file.path("output", "table"), recursive = T)

# Create the map
if (!file.exists(file.path("output/figure/map.pdf"))) {
  message("Plotting the map")
  system(paste("Rscript", file.path("code", "PlotDataMap.R")))
}

# Load libraries
message("loading the libraries")
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
  library(rlang)
  library(patchwork)
  library(tidygraph)
  library(ggraph)
  library(gridExtra)
  library(grid)
  #library(data.table)
  library(spaa)
  library(broom)
  library(ggrepel)
})
# Set the theme for plots
theme_set(
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.line = element_blank(),
      strip.background = element_blank()
    )
)
# Source the helper functions
source(file.path("code", "HelpMyAnalyses.R"))
# Import the data
message("importing the datasets")
timeseries_null <- readAndArrange(
  file = file.path("data", "analyses", "timeseries_null.csv")
) |>
  filter(year(sample_week) > 2007, isoweek(sample_week) %in% 2:51)
timeseries_fluxes <- readAndArrange(
  file = file.path("data", "analyses", "timeseries_fluxes.csv")
) |>
  filter(year(sample_week) > 2007, isoweek(sample_week) %in% 2:51)
timeseries_nodes <- readAndArrange(
  file = file.path("data", "analyses", "timeseries_nodes.csv")
) |>
  filter(year(sample_week) > 2007, isoweek(sample_week) %in% 2:51)
isoweek_fluxes <- readAndArrange(
  file = file.path("data", "analyses", "isoweek_fluxes.csv"),
  arrange = FALSE
) |>
  filter(iso_week %in% 2:51)
annual_fluxes <- readAndArrange(
  file = file.path("data", "analyses", "annual_fluxes.csv"),
  arrange = FALSE
)
station_fluxes <- readAndArrange(
  file = file.path("data", "analyses", "station_fluxes.csv"),
  arrange = FALSE
)
node_data <- readAndArrange(
  file = file.path("data", "raw", "node_data.csv"),
  arrange = FALSE
)
weekly_biomasses <- readAndArrange(
  file = file.path("data", "processed", "interpolation", "weekly_biomasses.csv")
) |>
  filter(year(sample_week) > 2007, isoweek(sample_week) %in% 2:51)
temperature <- readAndArrange(
  file = file.path("data", "processed", "interpolation", "temperature.csv")
)
# Set the color scheme for plotting
color_mapping <- setNames(node_data$color, node_data$node_name)

# Fig 1f ----
message("Generating Fig. 1f")
# Foodweb at the beggining of the timeseries and at the end
foodweb_graph <- timeseries_fluxes |>
  mutate(
    when = ifelse(
      sample_week == min(sample_week),
      "start",
      ifelse(sample_week == max(sample_week), "end", "none")
    )
  ) |>
  filter(when != "none") |>
  group_by(predator, prey, station, when) |>
  summarise(flux = mean(mean, na.rm = TRUE), .groups = "drop") |>
  group_by(
    TL = ifelse(
      predator %in% c("Clupea", "Sprattus", "Gasterosteus"),
      "fish",
      "zooplankton"
    ),
    when
  ) |>
  mutate(flux = (flux / sum(flux)) * 100) |>
  ungroup() |>
  as_tbl_graph() |>
  activate(nodes) |>
  left_join(rename(node_data, name = node_name), by = join_by(name)) |>
  activate(edges)
node_names <- foodweb_graph |> activate(nodes) |> pull(name)

# Extract edges and map names
edges_with_names <- foodweb_graph |>
  activate(edges) |>
  as_tibble() |>
  mutate(
    from_name = node_names[from],
    to_name = node_names[to]
  )

foodweb_graph <- foodweb_graph |>
  activate(edges) |>
  left_join(edges_with_names, by = join_by(from, to, station, when, flux, TL))

Fig1f <- ggraph(
  graph = foodweb_graph,
  layout = "manual",
  x = horizontal_position,
  y = trophic_level
) +
  # Add the link between all species relative to their contribution to the total fluxes between TL
  geom_edge_link(mapping = aes(edge_width = flux, col = to_name)) + #, alpha = flux)) +
  scale_edge_color_manual(values = color_mapping) +
  scale_edge_width(
    range = c(.5, 5),
    limits = c(0, 25),
    name = "Contribution fluxes [%]"
  ) +
  geom_point(
    data = node_data |> mutate(label = 1:24),
    mapping = aes(
      x = horizontal_position,
      y = trophic_level
    ),
    fill = "white",
    shape = 21,
    size = 6.5
  ) +
  geom_text(
    data = node_data |> mutate(label = 1:24),
    mapping = aes(
      x = horizontal_position,
      y = trophic_level,
      label = label
    )
  ) +
  scale_size_continuous(guide = "none") +
  facet_wrap(. ~ when) +
  coord_cartesian(xlim = c(0.8, 13.2), ylim = c(0.8, 3.2)) +
  scale_y_continuous(
    breaks = 1:3,
    labels = c("Phytoplankton", "Zooplankton", "Fish")
  ) +
  # Fix the theme
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL)
ggsave(
  plot = Fig1f,
  filename = file.path("output", "figure", "fig1f.pdf"),
  height = 3.5,
  width = 7,
  dpi = 500
)
# Average foodweb over the productive season ----
foodweb_avg <- timeseries_fluxes |>
  filter(
    station == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(predator, prey, station, season, year) |>
  summarise(flux = mean(mean, na.rm = TRUE), .groups = "drop_last") |>
  summarise(flux_avg = mean(flux, na.rm = TRUE), .groups = "drop") |>
  group_by(
    TL = ifelse(
      predator %in% c("Clupea", "Sprattus", "Gasterosteus"),
      "fish",
      "zooplankton"
    ),
    season
  ) |>
  mutate(
    rel_flux = (flux_avg / sum(flux_avg)) * 100,
    season = factor(season, levels = c("Spring", "Summer", "Fall"))
  ) |>
  ungroup() |>
  as_tbl_graph() |>
  activate(nodes) |>
  left_join(rename(node_data, name = node_name), by = join_by(name)) |>
  activate(edges)
node_names <- foodweb_avg |> activate(nodes) |> pull(name)

# Extract edges and map names
edges_with_names <- foodweb_avg |>
  activate(edges) |>
  as_tibble() |>
  mutate(
    from_name = node_names[from],
    to_name = node_names[to]
  )

foodweb_graph_avg <- foodweb_avg |>
  activate(edges) |>
  left_join(
    edges_with_names,
    by = join_by(from, to, station, season, flux_avg, rel_flux, TL)
  )

Fig3a <- ggraph(
  graph = foodweb_graph_avg,
  layout = "manual",
  x = horizontal_position,
  y = trophic_level
) +
  # Add the link between all species relative to their contribution to the total fluxes between TL
  geom_edge_link(mapping = aes(edge_width = rel_flux, col = to_name)) + #, alpha = flux)) +
  scale_edge_color_manual(values = color_mapping) +
  scale_edge_width(
    range = c(.2, 3),
    limits = c(0, 25),
    name = "Contribution fluxes [%]"
  ) +
  geom_point(
    data = node_data |> mutate(label = 1:24),
    mapping = aes(
      x = horizontal_position,
      y = trophic_level
    ),
    fill = "white",
    shape = 21,
    size = 6.5
  ) +
  geom_text(
    data = node_data |> mutate(label = 1:24),
    mapping = aes(
      x = horizontal_position,
      y = trophic_level,
      label = label
    )
  ) +
  scale_size_continuous(guide = "none") +
  coord_cartesian(xlim = c(0.8, 13.2), ylim = c(0.8, 3.2)) +
  scale_y_continuous(
    breaks = 1:3,
    labels = c("Phytoplankton", "Zooplankton", "Fish")
  ) +
  facet_grid(season ~ .) +
  # Fix the theme
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL, y = NULL)
ggsave(
  plot = Fig3a,
  filename = file.path("output", "figure", "fig3a.pdf"),
  height = 7.5,
  width = 4,
  dpi = 500
)

# Fig 2 ----
message("Generating Fig. 2")


# Seasonal biomasses
plot_data <- weekly_biomasses |>
  left_join(node_data, by = join_by(node_name)) |>
  filter(station_name == "BY31 LANDSORTSDJ", type != "fish") |>
  mutate(
    node_name = factor(node_name, levels = node_data$node_name),
    type = fct_rev(factor(
      type,
      levels = c("fish", "zooplankton", "phytoplankton")
    ))
  ) |>
  group_by(type, node_name, iso_week = isoweek(sample_week)) |>
  summarise(
    y = mean(biomass, na.rm = T),
    ymin = pmax(0, y - sd(biomass, na.rm = T)),
    ymax = y + sd(biomass, na.rm = T),
    .groups = "drop"
  )


# Add empty rows for 'fish'
fish_placeholder <- tibble(
  type = factor("fish", levels = c("fish", "zooplankton", "phytoplankton")),
  node_name = "Clupea",
  iso_week = 2,
  y = NA_real_,
  ymin = 0,
  ymax = 0
)

fig2.d <- weekly_biomasses |>
  left_join(node_data, by = join_by(node_name)) |>
  filter(station_name == "BY31 LANDSORTSDJ", type == "fish") |>
  mutate(
    node_name = factor(node_name, levels = node_data$node_name),
    type = fct_rev(factor(
      type,
      levels = c("fish", "zooplankton", "phytoplankton")
    ))
  ) |>
  group_by(type, year = year(sample_week), node_name) |>
  summarise(
    y = mean(biomass, na.rm = T),
    .groups = "drop_last"
  ) |>
  mutate(rel_biomass = y / sum(y)) |>
  ungroup() |>
  ggplot(aes(
    x = year,
    y = rel_biomass,
    fill = node_name
  )) +

  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_mapping) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "relative biomass", x = NULL)

Fig2.1 <- plot_data |>
  ggplot(aes(
    x = iso_week,
    y = y,
    ymin = ymin,
    ymax = ymax,
    col = node_name,
    fill = node_name
  )) +
  geom_ribbon(alpha = .1, colour = NA) +
  geom_line(linewidth = 1, na.rm = T) +
  facet_grid(type ~ ., scales = "free_y") +
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 20, 3),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Biomass [g/m2]")

Fig2.2 <- temperature |>
  filter(
    station_name == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  # add_season() |>
  #  filter(season != "Winter") |>
  group_by(iso_week = isoweek(sample_week)) |>
  summarise(
    y = mean(temperature, na.rm = T),
    ymin = y - sd(temperature, na.rm = T),
    ymax = y + sd(temperature, na.rm = T),
    type = "temperature",
    .groups = "drop"
  ) |>
  ggplot(aes(
    x = iso_week,
    y = y,
    ymin = ymin,
    ymax = ymax
  )) +
  geom_ribbon(alpha = .1, colour = NA) +
  geom_line(linewidth = 1) +
  facet_grid(type ~ .) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 20, 3),
    expand = c(0, 0),
    limits = c(0, 11)
  ) +
  labs(x = NULL, y = "Temperature")

Fig2.a <- Fig2.2 /
  Fig2.1 /
  fig2.d +
  plot_layout(heights = c(1, 2, 1), guides = "collect", axes = "collect") &
  guides(colour = guide_legend(ncol = 1))
bg_df <- tibble(
  year = 2008:2023,
  idx = seq_along(year)
) |>
  filter(idx %% 2 == 0)

#dodge_width <- 0.6
#fig2_colors <- c(
#  color_mapping,
#  Cladocerans = "#cfcecc",
#  Copepods = "#a84332",
#  Rotifers = "#8E518D",
#  Cyanobacteria = "#1c203b",
#  Diatoms = "#ede493",
#  Dinoflagellates = "#BD4A46",
#  Other = "#668704"
#)
# Prepare data and compute dodge positions
Fig2.3 <-
  weekly_biomasses |>
  left_join(node_data, by = join_by(node_name)) |>
  mutate(
    year = year(sample_week),
    node_name = factor(node_name, levels = node_data$node_name),
    type = factor(type, levels = c("phytoplankton", "zooplankton", "fish"))
  ) |>
  add_season() |>
  filter(
    station_name == "BY31 LANDSORTSDJ",
    season != "Winter"
  ) |>
  mutate(
    taxa = case_when(
      node_name %in% c("Bacillariales", "Chaetocerotales", "Thalassiosirales") ~
        "Diatoms",
      node_name == "Peridiniales" ~ "Dinoflagellates",
      node_name %in%
        c("Aphanizomenonaceae", "Nodularia", "Cyanobiaceae", "Pseudanabaena") ~
        "Cyanobacteria",
      node_name %in%
        c("Acartia", "Temora", "Eurytemora", "Centropages", "Pseudocalanus") ~
        "Copepods",
      node_name %in% c("Bosmina", "Evadne") ~ "Cladocerans",
      node_name == "Synchaeta" ~ "Rotifers",
      node_name %in% c("Sprattus", "Clupea", "Gasterosteus") ~ node_name,
      TRUE ~ "Other"
    )
  ) |>
  group_by(type, year, sample_week) |>
  summarise(biomass = sum(biomass), .groups = "drop_last") |>
  summarise(biomass = mean(biomass, na.rm = TRUE), .groups = "drop_last") |>
  mutate(
    z = (biomass - mean(biomass)) / sd(biomass),
    avg = mean(biomass)
  ) |>
  ungroup() |>
  # Manually compute dodge positions per year and type

  # group_by(type, year) |>
  #  mutate(
  #    n = n(),
  #    idx = row_number(),
  #    offset = ((idx - 1) / (n - 1)) - 0.5, # center around 0
  #    x_dodge = year + offset * dodge_width
  #  ) |>
  #  ungroup() |>

  # Plot
  ggplot(aes(x = year, y = biomass)) +
  # Background rectangles
  geom_rect(
    data = bg_df,
    aes(xmin = year - 0.5, xmax = year + 0.5, ymin = -Inf, ymax = Inf),
    fill = "gray80",
    inherit.aes = FALSE
  ) +
  geom_hline(mapping = aes(yintercept = avg), color = "black") +
  # Lollipop segments (manually dodged)
  geom_segment(
    aes(x = year, xend = year, y = avg, yend = biomass),
    linewidth = 0.4
  ) +
  # Lollipop heads
  geom_point(
    shape = 21,
    size = 3.5,
    color = "black",
    fill = "white"
  ) +
  scale_x_continuous(
    breaks = seq(2008, 2022, 2),
    expand = c(0, 0),
    limits = c(2007.5, 2023.5)
  ) +
  facet_grid(type ~ ., scales = "free_y") +
  labs(x = NULL, y = "Biomass") # +


Fig2.4 <- temperature |>
  filter(
    station_name == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(year = year(sample_week)) |>
  summarise(temp_avg = mean(temperature, na.rm = T), .groups = "drop") |>
  #  group_by(season) |>
  mutate(
    type = "temperature",
    z = (temp_avg - mean(temp_avg)) / sd(temp_avg),
    avg = mean(temp_avg),
    col = ifelse(z < 0, "neg", "pos")
  ) |>
  ggplot(aes(x = year, y = temp_avg)) +

  geom_rect(
    data = bg_df,
    aes(
      xmin = year - 0.5,
      xmax = year + 0.5,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "gray80",
    inherit.aes = FALSE
  ) +
  geom_hline(mapping = aes(yintercept = avg), color = "black") +
  geom_segment(
    aes(x = year, xend = year, y = avg, yend = temp_avg),
    linewidth = 0.4
  ) +

  # Lollipop heads
  geom_point(
    shape = 21,
    size = 3,
    color = "black",
    fill = "white"
  ) +
  #scale_fill_manual(values = c("#848FA5", "#C14953"), guide = "none") +
  scale_x_continuous(
    breaks = seq(2008, 2022, 2),
    expand = c(0, 0),
    limits = c(2007.5, 2023.5)
  ) +

  labs(x = NULL, y = "Temperature") +
  facet_grid(type ~ .)
Fig2.b <- Fig2.4 /
  Fig2.3 +
  plot_layout(heights = c(1, 3), guides = "collect", axes = "collect")
ggsave(
  plot = Fig2.a,
  filename = file.path("output", "figure", "fig2a.pdf"),
  width = 8,
  height = 8
)
ggsave(
  plot = Fig2.b,
  filename = file.path("output", "figure", "fig2b.pdf"),
  width = 5,
  height = 8
)
# Add the stats ....
df_test_temperature <- temperature |>
  filter(
    station_name == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(year = year(sample_week)) |>
  summarise(temp_avg = mean(temperature, na.rm = T), .groups = "drop") |>
  #  group_by(season) |>
  mutate(
    type = "temperature",
    z = (temp_avg - mean(temp_avg)) / sd(temp_avg),
    param = temp_avg
  )
df_test_biomass <- weekly_biomasses |>
  left_join(node_data, by = join_by(node_name)) |>
  mutate(
    year = year(sample_week),
    node_name = factor(node_name, levels = node_data$node_name),
    type = factor(type, levels = c("phytoplankton", "zooplankton", "fish"))
  ) |>
  add_season() |>
  filter(
    station_name == "BY31 LANDSORTSDJ",
    season != "Winter"
  ) |>
  group_by(type, year, sample_week) |>
  summarise(biomass = sum(biomass), .groups = "drop_last") |>
  summarise(biomass = mean(biomass, na.rm = TRUE), .groups = "drop_last") |>
  mutate(
    z = (biomass - mean(biomass)) / sd(biomass)
  ) |>
  ungroup() |>
  mutate(param = biomass)

run_model_workflow(
  df = bind_rows(df_test_biomass, df_test_temperature),
  group_var = type,
  formula_expr = "param ~ year",
  filename_suffix = "_vs_year.pdf",
  plot_title = "vs year",
  slope = "year"
) |>
  select(type, term, estimate, r.squared, p.value) |>
  filter(term == "slope") |>
  ungroup() |>
  mutate(p_adjusted = p.adjust(p.value, method = "fdr")) |>
  write.csv(file = file.path("output", "table", "trends.csv"))
# Fig 2.5 ----
message("Generating figure 2.5")
fluxes <- timeseries_fluxes |>
  filter(
    station == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  left_join(node_data, by = join_by(predator == node_name)) |>
  group_by(
    year = year(sample_week),
    isoweek = isoweek(sample_week),
    sample_week,
    prey,
    type,
    predator
  ) |>
  summarise(flux = mean(mean, na.rm = T), .groups = "drop_last") |>
  summarise(flux = sum(flux), .groups = "drop") |>
  group_by(type, year, isoweek) |>
  summarise(flux = sum(flux), .groups = "drop") |>
  mutate(prey = ifelse(type == "fish", "zooplankton", "phytoplankton"))
biomass <- weekly_biomasses |>
  filter(
    station_name == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  left_join(node_data, by = join_by(node_name)) |>
  filter(type != "fish") |>
  group_by(
    year = year(sample_week),
    isoweek = isoweek(sample_week),
    sample_week,
    type
  ) |>
  summarise(biomass = sum(biomass, na.rm = T), .groups = "drop") |>
  rename("prey" = type)
predation_pressure_df <- biomass |>
  left_join(fluxes, by = join_by(year, isoweek, prey)) |>
  mutate(
    predation_pressure = flux / biomass,
    type = prey,
    predation_pressure = ifelse(
      type == "zooplankton",
      predation_pressure * 5,
      predation_pressure
    )
  ) |>
  group_by(isoweek, type, parameter = "predation_pressure") |>
  summarise(
    y = mean(predation_pressure),
    ymin = pmax(0, y - sd(predation_pressure)),
    ymax = y + sd(predation_pressure),
    .groups = "drop"
  )

consumption_df <-
  fluxes |>
  mutate(
    type = ifelse(type == "zooplankton", "phytoplankton", "zooplankton"),
    flux = ifelse(type == "zooplankton", flux * 10, flux)
  ) |>
  group_by(isoweek, type, parameter = "consumed") |>
  summarise(
    y = mean(flux),
    ymin = pmax(0, y - sd(flux)),
    ymax = y + sd(flux),
    .groups = "drop"
  )
ratio_df <- fluxes |>
  group_by(isoweek, type, year) |>
  summarise(flux_tot = sum(flux), .groups = "drop_last") |>
  pivot_wider(names_from = type, values_from = flux_tot) |>
  mutate(ratio = fish / zooplankton) |>
  group_by(isoweek, type = "ratio", parameter = "ratio") |>
  summarise(
    y = mean(ratio),
    ymin = y - sd(ratio),
    ymax = y + sd(ratio),
    .groups = "drop"
  )

biomass_df <- biomass |>
  mutate(type = prey) |>
  group_by(isoweek, type, parameter = "biomass") |>
  summarise(
    y = mean(biomass),
    ymin = pmax(0, y - sd(biomass)),
    ymax = y + sd(biomass),
    .groups = "drop"
  )

fig2.5 <- bind_rows(
  consumption_df,
  ratio_df,
  predation_pressure_df #,
  #biomass_df
) |>
  ggplot(aes(
    x = isoweek,
    y = y,
    ymin = ymin,
    ymax = ymax,
    fill = type,
    col = type
  )) +
  geom_ribbon(alpha = .2, col = NA) +
  geom_line() +
  facet_grid(parameter ~ ., scales = "free_y") +
  scale_fill_manual(
    values = c(
      "phytoplankton" = "#006837",
      "zooplankton" = "#fe9929",
      "ratio" = "black"
    )
  ) +
  scale_color_manual(
    values = c(
      "phytoplankton" = "#006837",
      "zooplankton" = "#fe9929",
      "ratio" = "black"
    )
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb, #c("Jan", "Mar", "May", "Jul", "Sep", "Nov"),
    expand = c(0, 0)
  ) +
  labs(x = NULL) +
  theme(legend.position = "bottom")


pred_pressure_ts <- biomass |>
  left_join(fluxes, by = join_by(year, isoweek, prey)) |>
  mutate(
    predation_pressure = flux / biomass,
    type = prey,
    predation_pressure = ifelse(
      type == "zooplankton",
      predation_pressure * 5,
      predation_pressure
    )
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(year, type, parameter = "predation_pressure") |>
  summarise(
    y = mean(predation_pressure),
    .groups = "drop"
  )

flux_ts <- timeseries_fluxes |>
  filter(isoweek(sample_week) %in% 2:51, year %in% 2008:2023) |>
  add_season() |>
  filter(season != "Winter") |>
  select(-c(upper, lower)) |>
  pivot_wider(names_from = prey, values_from = mean, values_fill = 0) |>
  pivot_longer(7:26, names_to = "prey", values_to = "flux") |>
  group_by(predator, prey) |>
  filter(sum(flux) > 0) |>
  ungroup() |>
  group_by(predator, prey, year, iso_week) |>
  summarise(weekly_flux = mean(flux), .groups = "drop_last") |>
  summarise(
    n_week = n_distinct(iso_week),
    n_days = n_week * 7,
    flux = sum(weekly_flux, na.rm = T) * 7,
    .groups = "drop"
  ) |>
  left_join(node_data, by = join_by(prey == node_name)) |>
  group_by(type, parameter = "consumed", year) |>
  summarise(y = sum(flux), .groups = "drop") |>
  mutate(
    y = ifelse(type == "zooplankton", y * 10, y)
  )

ratio_ts <- timeseries_fluxes |>
  filter(
    station == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  left_join(node_data, by = join_by(predator == node_name)) |>
  group_by(
    year = year(sample_week),
    isoweek = isoweek(sample_week),
    sample_week,
    prey,
    type,
    predator
  ) |>
  summarise(flux = mean(mean, na.rm = T), .groups = "drop_last") |>
  summarise(flux = sum(flux), .groups = "drop") |>
  group_by(type, year, isoweek, sample_week) |>
  summarise(flux = sum(flux), .groups = "drop") |>
  mutate(prey = ifelse(type == "fish", "zooplankton", "phytoplankton")) |>
  group_by(isoweek, sample_week, type, year) |>
  summarise(flux_tot = sum(flux), .groups = "drop_last") |>
  pivot_wider(names_from = type, values_from = flux_tot) |>
  mutate(ratio = fish / zooplankton) |>
  ungroup() |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(year, type = "ratio", parameter = "ratio") |>
  summarise(
    y = mean(ratio),
    .groups = "drop"
  )
fig2.6 <- bind_rows(ratio_ts, pred_pressure_ts, flux_ts) |>
  ggplot(aes(x = year, y = y, fill = type, col = type)) +
  geom_line() +
  geom_point(shape = 21, col = "black", size = 2) +
  scale_fill_manual(
    values = c(
      "phytoplankton" = "#006837",
      "zooplankton" = "#fe9929",
      "ratio" = "black"
    )
  ) +
  scale_color_manual(
    values = c(
      "phytoplankton" = "#006837",
      "zooplankton" = "#fe9929",
      "ratio" = "black"
    )
  ) +
  theme(legend.position = "bottom") +
  facet_grid(parameter ~ ., scales = "free")
new_fig2.5 <- fig2.6 + fig2.5
ggsave(
  plot = new_fig2.5,
  filename = file.path("output", "figure", "fig2.5.pdf"),
  width = 8,
  height = 6.5
)
# Fig 3 ----
message("Generating Fig. 3")
## Relative flux by week ----
tot_fluxes <- timeseries_fluxes |>
  group_by(
    sample_week,
    predator = factor(predator, levels = node_data$node_name),
    station
  ) |>
  summarise(tot_flux = sum(mean, na.rm = T), .groups = "drop") |>
  filter(isoweek(sample_week) %in% 2:51) |>
  group_by(iso_week = isoweek(sample_week), predator) |>
  summarise(
    consum = mean(tot_flux),
    .groups = "drop"
  )

# Average contribution of prey per week:
rel_contribution_fluxes <-
  timeseries_fluxes |>
  group_by(sample_week, predator) |>
  reframe(
    sample_week = sample_week,
    station = station,
    prey = factor(prey, levels = node_data$node_name),
    predator = factor(predator, levels = node_data$node_name),
    type = ifelse(
      predator %in% c("Clupea", "Gasterosteus", "Sprattus"),
      "fish",
      "zooplankton"
    ),
    rel_flux = mean / sum(mean, na.rm = T)
  ) |>
  ungroup()

seasonal_fluxes <- rel_contribution_fluxes |>
  pivot_wider(names_from = prey, values_from = rel_flux, values_fill = 0) |>
  pivot_longer(where(is.numeric), names_to = "prey", values_to = "rel_flux") |>
  filter(isoweek(sample_week) %in% 2:51) |>
  group_by(iso_week = isoweek(sample_week), predator, prey) |>
  summarise(
    y = mean(rel_flux, na.rm = T),
    .groups = "drop"
  ) |>
  mutate(
    predator = factor(predator, levels = node_data$node_name),
    prey = factor(prey, levels = node_data$node_name)
  ) |>
  left_join(tot_fluxes, by = join_by(iso_week, predator)) |>
  mutate(
    contribution = consum * y,
    type = ifelse(
      predator %in% c("Clupea", "Sprattus", "Gasterosteus"),
      "fish",
      "zooplankton"
    ),
    type = factor(type, levels = c("zooplankton", "fish"))
  )
# Average consumption over the entire time series during the productive season
flux_productive_season <- timeseries_fluxes |>
  filter(isoweek(sample_week) %in% 2:51, year %in% 2008:2023) |>
  add_season() |>
  filter(season != "Winter") |>
  select(-c(upper, lower)) |>
  pivot_wider(names_from = prey, values_from = mean, values_fill = 0) |>
  pivot_longer(7:26, names_to = "prey", values_to = "flux") |>
  group_by(predator, prey) |>
  filter(sum(flux) > 0) |>
  ungroup() |>
  group_by(predator, prey, iso_week) |>
  summarise(weekly_flux = mean(flux), .groups = "drop_last") |>
  summarise(
    n_week = n_distinct(iso_week),
    n_days = n_week * 7,
    flux = sum(weekly_flux, na.rm = T) * 7,
    .groups = "drop"
  )
order_predator <- flux_productive_season |>
  group_by(predator) |>
  summarise(tot = sum(flux), .groups = "drop") |>
  arrange(tot) |>
  pull(predator)
order_prey <- flux_productive_season |>
  group_by(prey) |>
  summarise(tot = sum(flux), .groups = "drop") |>
  arrange(tot) |>
  pull(prey)

fig3.2 <-
  flux_productive_season |>
  mutate(
    predator = fct_rev(factor(predator, levels = order_predator)),
    prey = factor(prey, levels = node_data$node_name),
    type = ifelse(
      predator %in% c("Gasterosteus", "Clupea", "Sprattus"),
      "2_fish",
      "1_zooplankton"
    ),
    flux = ifelse(type == "2_fish", flux * 10, flux)
  ) |>
  ggplot(aes(y = flux, x = predator, fill = prey)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +
  facet_grid(. ~ type, scales = "free", space = "free") +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 25),
    breaks = seq(0, 24, 5),
    sec.axis = sec_axis(transform = ~ . / 10)
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(x = NULL, y = "Total consumption from spring to fall")
fig3.3 <- flux_productive_season |>
  mutate(
    predator = factor(predator, levels = node_data$node_name),
    prey = fct_rev(factor(prey, levels = order_prey)),
    type = ifelse(
      predator %in% c("Gasterosteus", "Clupea", "Sprattus"),
      "2_fish",
      "1_zooplankton"
    ),
    flux = ifelse(type == "2_fish", flux * 10, flux)
  ) |>
  ggplot(aes(y = flux, fill = predator, x = prey)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(
    expand = c(0, 0),
    #limits = c(0, 22),
    #breaks = seq(0, 22, 5)
    sec.axis = sec_axis(transform = ~ . / 10)
  ) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  facet_grid(. ~ type, scales = "free", space = "free") +
  labs(x = NULL, y = "Total outgoing fluxes from spring to fall")


fig3.1 <-
  seasonal_fluxes |>
  mutate(predator = fct_rev(factor(predator, levels = order_predator))) |>
  ggplot(aes(x = iso_week, y = contribution, fill = prey)) +

  geom_area(stat = "identity") +
  facet_grid(predator ~ ., scales = "free_y") +
  scale_fill_manual(values = color_mapping) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb, #c("Jan", "Mar", "May", "Jul", "Sep", "Nov"),
    expand = c(0, 0)
  ) +
  labs(x = NULL, y = "Consumption \n [kJ/day/m2]") +
  theme(legend.position = "right")
# Consumption annomalies
fig3.4 <- timeseries_fluxes |>
  filter(isoweek(sample_week) %in% 2:51, year %in% 2008:2023) |>
  add_season() |>
  filter(season != "Winter") |>
  select(-c(upper, lower)) |>
  pivot_wider(names_from = prey, values_from = mean, values_fill = 0) |>
  pivot_longer(7:26, names_to = "prey", values_to = "flux") |>
  group_by(predator, prey) |>
  filter(sum(flux) > 0) |>
  ungroup() |>
  group_by(predator, year, iso_week) |>
  mutate(flux = sum(flux)) |>
  summarise(weekly_flux = mean(flux), .groups = "drop_last") |>
  summarise(
    flux = sum(weekly_flux, na.rm = T),
    .groups = "drop"
  ) |>
  group_by(predator) |>
  mutate(
    anomalie = (flux - mean(flux)) / (sd(flux)),
    predator = factor(predator, levels = order_predator)
  ) |>
  ungroup() |>
  ggplot(aes(x = year, xend = year, y = anomalie, yend = 0)) +
  geom_hline(yintercept = 0) +

  geom_segment() +
  geom_point(shape = 21, fill = "white", size = 2) +
  facet_grid(predator ~ .) +
  scale_x_continuous(breaks = seq(2008, 2023, 2)) +
  scale_y_continuous(limits = c(-3, 3))


fig3 <- fig3.1 +
  (fig3.2 /
    fig3.3) +
  plot_layout(guides = "collect", widths = c(1, 2))

ggsave(
  plot = fig3,
  filename = file.path("output", "figure", "fig3.pdf"),
  width = 15,
  height = 10
)

# Fig S13 ----
message("Generating Sup. Fig. S13")
df_pp <- timeseries_fluxes |>
  mutate(
    predator = factor(predator, levels = node_data$node_name),
    prey = factor(prey, levels = node_data$node_name),
    type = if_else(
      predator %in% c("Clupea", "Gasterosteus", "Sprattus"),
      "fish",
      "zooplankton"
    )
  ) |>
  group_by(sample_week, station, predator, prey, type) |>
  summarise(flux = sum(mean, na.rm = TRUE), .groups = "drop") |>
  filter(type == "fish") |>

  # Join with zooplankton → phytoplankton contribution
  left_join(
    rel_contribution_fluxes |>
      filter(type == "zooplankton", rel_flux > 0) |>
      pivot_wider(names_from = prey, values_from = rel_flux, values_fill = 0) |>
      rename(prey = predator),
    by = join_by(sample_week, prey, station)
  ) |>

  # Scale fish flux by zooplankton contributions to phytoplankton
  mutate(across(where(is.numeric) & !any_of("flux"), ~ .x * flux)) |>
  select(-flux) |>
  pivot_longer(
    where(is.numeric),
    names_to = "primary_producer",
    values_to = "scaled_rel_flux"
  ) |>

  # Aggregate to fish ← phytoplankton step
  group_by(sample_week, predator, primary_producer, station) |>
  summarise(
    rel_primary_production = mean(scaled_rel_flux, na.rm = TRUE),
    .groups = "drop"
  ) |>

  # Normalize within fish predator groups
  group_by(sample_week, predator, station) |>
  mutate(
    rel_primary_production = rel_primary_production /
      sum(rel_primary_production)
  ) |>
  ungroup()


SupFigS13.1 <- df_pp |>
  mutate(
    iso_week = isoweek(sample_week),
    predator = factor(predator, levels = node_data$node_name),
    primary_producer = factor(primary_producer, levels = node_data$node_name)
  ) |>
  filter(iso_week %in% 2:51) |>
  group_by(iso_week, predator, primary_producer) |>
  summarise(y = mean(rel_primary_production, na.rm = TRUE), .groups = "drop") |>

  left_join(tot_fluxes, by = join_by(iso_week, predator)) |>
  mutate(
    contribution = consum * y,
    type = if_else(
      predator %in% c("Clupea", "Sprattus", "Gasterosteus"),
      "fish",
      "zooplankton"
    ),
    type = factor(type, levels = c("zooplankton", "fish"))
  ) |>

  ggplot(aes(x = iso_week, y = contribution, fill = primary_producer)) +
  geom_area() +
  facet_grid(predator ~ type, scales = "fixed") +
  scale_fill_manual(values = color_mapping) +
  scale_x_continuous(
    breaks = seq(1, 52, 4.35),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  scale_y_continuous(breaks = seq(0, 0.02, 0.005)) +
  labs(x = NULL, y = "Consumption \n [kJ/day/m2]")

phyto_to_fish <- df_pp |>
  mutate(
    iso_week = isoweek(sample_week),
    year = year(sample_week)
  ) |>
  filter(iso_week %in% 2:51, year %in% 2008:2023) |>
  add_season() |>
  filter(season != "Winter") |>

  group_by(predator, prey = primary_producer, iso_week) |>
  summarise(
    y = mean(rel_primary_production, na.rm = TRUE),
    .groups = "drop_last"
  ) |>
  summarise(
    n_week = n_distinct(iso_week),
    n_days = n_week * 7,
    rel_primary_production = mean(y, na.rm = TRUE),
    .groups = "drop"
  ) |>

  group_by(predator) |>
  mutate(
    rel_primary_production = rel_primary_production /
      sum(rel_primary_production)
  ) |>
  ungroup() |>

  left_join(
    flux_productive_season |>
      group_by(predator) |>
      summarise(tot_flux = sum(flux), .groups = "drop"),
    by = join_by(predator)
  )

SupFigS13.2 <- phyto_to_fish |>
  mutate(
    phytoplankton_contribution = rel_primary_production * tot_flux,
    predator = fct_rev(factor(predator, levels = order_predator)),
    prey = factor(prey, levels = node_data$node_name)
  ) |>

  ggplot(aes(y = phytoplankton_contribution, x = predator, fill = prey)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping, guide = "none") +
  scale_y_continuous(expand = c(0, 0), breaks = 0:2) +
  #facet_grid(predator ~ ., scales = "free_y") +
  theme(legend.position = "right") +
  labs(x = "Total consumption from spring to fall", y = NULL)
order_phyto <-
  phyto_to_fish |>
  mutate(
    phytoplankton_contribution = rel_primary_production * tot_flux,
  ) |>
  group_by(prey) |>
  summarise(tot = sum(phytoplankton_contribution), .groups = "drop") |>
  arrange(tot) |>
  pull(prey)
SupFigS13.3 <- phyto_to_fish |>
  mutate(
    phytoplankton_contribution = rel_primary_production * tot_flux,
    predator = fct_rev(factor(predator, levels = node_data$node_name)),
    prey = fct_rev(factor(prey, levels = order_phyto))
  ) |>

  ggplot(aes(y = phytoplankton_contribution, fill = predator, x = prey)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_mapping, guide = "none") +
  scale_y_continuous(expand = c(0, 0), breaks = 0:2) +
  #facet_grid(predator ~ ., scales = "free_y") +
  theme(legend.position = "right") +
  labs(x = "Total outgoing fluxes from spring to fall", y = NULL)

SupFigS13 <-
  SupFigS13.1 +
  (SupFigS13.2 /
    SupFigS13.3 &
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))) +
  plot_layout(guides = "collect", width = c(1, 2))
ggsave(
  plot = SupFigS13,
  filename = file.path("output", "figure", "SupFigS13.pdf"),
  width = 10,
  height = 6
)

# Ordination ----
message("Ordination in progress...")
# For phytoplankton
mat_pp <- df_pp |>
  pivot_wider(
    names_from = primary_producer,
    values_from = rel_primary_production,
    values_fill = 0
  ) |>
  na.omit()
bray_dist_pp <- vegdist(select(mat_pp, where(is.numeric)), method = "bray")
pcoa_result_pp <- pcoa(bray_dist_pp)
site_scores_pp <- pcoa_result_pp$vectors
envfit_pp <- as.data.frame(scores(
  envfit(site_scores_pp, select(mat_pp, where(is.numeric)), permutations = 999),
  display = c("vectors")
)) |>
  rownames_to_column(var = "prey") |>
  rename("Axis1" = Axis.1, "Axis2" = Axis.2)
eig_pp <- pcoa_result_pp$values$Relative_eig
pcoa_df_pp <- process_pcoa(pcoa_result_pp, metadata = mat_pp)
ggsave(
  plot = PCOA_plot(pcoa_df_pp, eig_pp),
  filename = file.path("output", "figure", "SupFigS15.pdf"),
  dpi = 500,
  width = 8,
  height = 6
)

# For zooplankton
mat_zp <-
  rel_contribution_fluxes |>
  filter(type == "fish", rel_flux > 0) |>
  select(sample_week, station, predator, prey, rel_flux) |>
  pivot_wider(names_from = prey, values_from = rel_flux, values_fill = 0)
bray_dist_zp <- vegdist(select(mat_zp, where(is.numeric)), method = "bray")
pcoa_result_zp <- pcoa(bray_dist_zp)
site_scores_zp <- pcoa_result_zp$vectors
envfit_zp <- as.data.frame(scores(
  envfit(site_scores_zp, select(mat_zp, where(is.numeric)), permutations = 999),
  display = c("vectors")
)) |>
  rownames_to_column(var = "prey") |>
  rename("Axis1" = Axis.1, "Axis2" = Axis.2)
eig_zp <- pcoa_result_zp$values$Relative_eig
pcoa_df_zp <- process_pcoa(pcoa_result_zp, metadata = mat_zp)
ggsave(
  plot = PCOA_plot(pcoa_df_zp, eig_zp),
  filename = file.path("output", "figure", "SupFigS14.pdf"),
  dpi = 500,
  width = 8,
  height = 6
)
# bind envfit with type
envfit <- bind_rows(
  envfit_zp |> mutate(type = "zooplankton"),
  envfit_pp |> mutate(type = "phytoplankton")
)
# Standardized flux
standardized_flux <- tot_fluxes |>
  group_by(predator) |>
  mutate(z = scale(consum)[, 1]) |>
  ungroup()

# Plot settings
shape_vals <- c(21, 22, 24)
fill_palette <- c("#F3C178", "#D8F1A0", "#065C58")

# Fig 4.1 ----
message("Generating Fig. 4")
Fig4.1 <- bind_pcoa(pcoa_df_pp, pcoa_df_zp) |>
  group_by(iso_week, predator, season, type) |>
  summarise(
    avg1 = mean(Axis1, na.rm = TRUE),
    avg2 = mean(Axis2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(standardized_flux, by = join_by(iso_week, predator)) |>
  filter(iso_week %in% 2:51) |>
  mutate(type = factor(type, levels = c("zooplankton", "phytoplankton"))) |>
  ggplot() +
  geom_segment(
    data = envfit |>
      mutate(type = factor(type, levels = c("zooplankton", "phytoplankton"))),
    aes(x = Axis1, y = Axis2, xend = 0, yend = 0),
    arrow = arrow(length = unit(2, 'mm'), ends = "first")
  ) +
  geom_text_repel(
    data = envfit |>
      mutate(type = factor(type, levels = c("zooplankton", "phytoplankton"))),
    aes(x = Axis1, y = Axis2, label = prey),
    size = 3,
    seed = 100
  ) +
  geom_path(aes(x = avg1, y = avg2, col = predator, linewidth = z)) +
  geom_point(
    aes(
      x = avg1,
      y = avg2,
      fill = iso_week,
      col = predator,
      group = predator,
      shape = predator,
      size = z
    )
  ) +
  facet_grid(. ~ type) +
  coord_fixed() +
  scale_shape_manual(values = shape_vals) +
  scale_color_manual(values = color_mapping, guide = "none") +
  scale_fill_gradientn(
    colors = fill_palette,
    values = scales::rescale(c(2, 22, 48)),
    limits = c(1, 52),
    breaks = c(12, 24, 36, 48),
    labels = c("Mar", "Jun", "Sep", "Dec")
  ) +
  scale_size_continuous(
    guide = "none",
    limits = c(-1.5, 1.6),
    range = c(1, 5)
  ) +
  scale_linewidth_continuous(
    guide = "none",
    limits = c(-1.5, 1.6),
    range = c(1, 5)
  ) +
  scale_x_continuous(limits = c(-0.9, 0.9)) +
  scale_y_continuous(limits = c(-0.9, 0.9)) +
  labs(x = "PCOA1", y = "PCOA2")


# -------
bind_pcoa(pcoa_df_pp, pcoa_df_zp) |>
  filter(season != "Winter") |>
  group_by(year = year(sample_week), predator, type) |>
  summarise(
    pcoa1 = mean(Axis1, na.rm = TRUE),
    pcoa2 = mean(Axis2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(4:5, names_to = "axis", values_to = "pcoa") |>
  ggplot(aes(x = year, y = pcoa, fill = predator, shape = predator)) +
  geom_point(size = 3) +
  facet_grid(type ~ axis) +
  scale_shape_manual(values = shape_vals)
#scale_fill_manual(values = color_mapping)
temp_data_standard <- temperature |>
  filter(
    station_name == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(year = year(sample_week)) |>
  summarise(temp_avg = mean(temperature, na.rm = T), .groups = "drop") |>
  #  group_by(season) |>
  mutate(
    z = (temp_avg - mean(temp_avg)) / sd(temp_avg),
    col = ifelse(z < 0, "neg", "pos")
  )
bind_pcoa(pcoa_df_pp, pcoa_df_zp) |>
  filter(season != "Winter", type == "zooplankton") |>
  group_by(year = year(sample_week), predator, type) |>
  summarise(
    pcoa1 = mean(Axis1, na.rm = TRUE),
    pcoa2 = mean(Axis2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(4:5, names_to = "axis", values_to = "pcoa") |>
  left_join(temp_data_standard, by = join_by(year)) |>
  ggplot(aes(x = z, y = pcoa, fill = predator, shape = predator)) +
  geom_point(size = 3) +
  facet_grid(type ~ axis) +
  scale_shape_manual(values = shape_vals) +
  geom_smooth(method = "lm", aes(col = predator), se = F) +
  labs(x = "temperature anomalies", y = NULL)
#Fig4 <- Fig4.1 /
#  Fig4.2 +
#  plot_layout(axis_titles = "collect", axes = "collect")
ggsave(
  plot = Fig4.1,
  filename = file.path("output", "figure", "fig4.pdf"),
  dpi = 500,
  width = 10,
  height = 5
)

# Table S3: Permanovas ----
# This step takes some time about 5 min so it is better to save the output to come back to it later
if (file.exists(file.path("output", "table", "permanova_zp.csv"))) {
  permanova_zp <- read_csv(
    file = file.path("output", "table", "permanova_zp.csv"),
    show_col_types = FALSE
  )
} else {
  message(
    "Running the permANOVA to test whether fish diet composition varies with time and fish species \n This takes about 5 min, you have plenty of time to get a coffee..."
  )
  set.seed(100)
  permanova_zp <- adonis2(
    formula = bray_dist_zp ~ predator + factor(iso_week),
    data = mutate(mat_zp, iso_week = isoweek(sample_week)),
    permutations = 999,
    #strata = mat_zp$month,
    parallel = 8,
    by = "margin"
  )
  permanova_zp |>
    as.data.frame() |>
    rownames_to_column("term") |>
    as_tibble() |>
    write_csv(file = file.path("output", "table", "permanova_zp.csv"))
}

if (file.exists(file.path("output", "table", "permanova_pp.csv"))) {
  permanova_pp <- read_csv(
    file = file.path("output", "table", "permanova_pp.csv"),
    show_col_types = FALSE
  )
} else {
  message(
    "Running the permANOVA to test whether primary producer composition varies with time and predator \n This takes about 5 min, you have plenty of time to get a coffee... \n If you already had your coffee, maybe this is a sign to force you to take a well deserved break..."
  )
  set.seed(100)
  permanova_pp <- adonis2(
    formula = bray_dist_pp ~ predator + factor(iso_week),
    data = mutate(mat_pp, iso_week = isoweek(sample_week)),
    permutations = 999,
    #strata = mat_pp$month,
    parallel = 8,
    by = "margin"
  )
  permanova_pp |>
    as.data.frame() |>
    rownames_to_column("term") |>
    as_tibble() |>
    write_csv(file = file.path("output", "table", "permanova_pp.csv"))
}

# Fig S16 -----
message("Generating Sup. Fig. S16")
inter_overlap_zp <-
  mat_zp |>
  filter(isoweek(sample_week) %in% 2:51) |>
  mutate(
    group_id = paste(year(sample_week), isoweek(sample_week), sep = "_")
  ) |>
  group_by(group_id) |>
  group_split() |>
  map_dfr(
    ~ interspecific_overlap(.x, 4:11)
  ) |>
  mutate(facet = "interspecific", type = "zp") |>
  filter(x != y)

inter_overlap_pp <- mat_pp |>
  filter(isoweek(sample_week) %in% 2:51) |>
  mutate(
    group_id = paste(year(sample_week), isoweek(sample_week), sep = "_")
  ) |>
  group_by(group_id) |>
  group_split() |>
  map_dfr(
    ~ interspecific_overlap(.x, 4:15)
  ) |>
  mutate(facet = "interspecific") |>
  filter(x != y)

SupFigS16.a <- inter_overlap_zp |>
  bind_rows(inter_overlap_pp |> mutate(type = "pp")) |>
  mutate(
    interaction = paste(x, y, sep = "_"),
    type = factor(type, levels = c("zp", "pp"))
  ) |>
  group_by(iso_week = isoweek(sample_week.x), interaction, facet, type) |>
  summarise(
    y = mean(value),
    ymin = pmax(0, y - sd(value)),
    ymax = pmin(1, y + sd(value)),
    .groups = "drop"
  ) |>
  ggplot(aes(
    x = iso_week,
    y = y,
    ymin = ymin,
    ymax = ymax,
    col = interaction,
    fill = interaction,
    group = interaction
  )) +
  geom_ribbon(alpha = .2, linewidth = 0) +
  geom_line(linewidth = 1) +
  labs(x = NULL, y = "Schoener's D") +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = c(
      "Clupea_Gasterosteus" = "#74a9cf",
      "Clupea_Sprattus" = "#ffffd4",
      "Gasterosteus_Sprattus" = "#b30000"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Clupea_Gasterosteus" = "#74a9cf",
      "Clupea_Sprattus" = "black",
      "Gasterosteus_Sprattus" = "#b30000"
    )
  ) +

  facet_grid(type ~ ., scales = "fixed")

SupFig16.b <- inter_overlap_zp |>
  bind_rows(inter_overlap_pp |> mutate(type = "pp")) |>
  mutate(
    interaction = paste(x, y, sep = "_"),
    type = factor(type, levels = c("zp", "pp"))
  ) |>
  mutate(year = year(sample_week.x), iso_week = isoweek(sample_week.x)) |>
  filter(iso_week %in% 11:48) |>
  group_by(year, interaction, type) |>
  summarise(
    y = mean(value),
    .groups = "drop"
  ) |>
  ggplot(aes(
    x = year,
    y = y,
    col = interaction,
    fill = interaction,
    group = interaction
  )) +
  geom_line(linewidth = 1) +
  geom_point(shape = 21, color = "black", size = 3) +
  labs(x = NULL, y = "Schoener's D") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(
    values = c(
      "Clupea_Gasterosteus" = "#74a9cf",
      "Clupea_Sprattus" = "black",
      "Gasterosteus_Sprattus" = "#b30000"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Clupea_Gasterosteus" = "#74a9cf",
      "Clupea_Sprattus" = "#ffffd4",
      "Gasterosteus_Sprattus" = "#b30000"
    )
  ) +
  facet_grid(type ~ ., scales = "fixed")
SupFigS16 <- SupFig16.b + SupFigS16.a + plot_layout(guides = "collect")
ggsave(
  plot = SupFigS16,
  filename = file.path("output", "figure", "SupFigS16.pdf"),
  width = 9,
  height = 4
)
# Fig S17-18 ----
message("Generating Sup. Fig. S17 and S18")
# Impact of forage ratios
# or how the forage ratios deviate from the null model
isoweek_null <- timeseries_null |>
  filter(year(sample_week) > 2007) |>
  group_by(predator, prey, iso_week = isoweek(sample_week), station) |>
  summarise(avg_fluxes = mean(flux, na.rm = T), .groups = "drop")

ggsave(
  filename = file.path("output", "figure", "SupFigS17.pdf"),
  plot = plot_flux_difference("phytoplankton"),
  height = 10,
  width = 8,
  dpi = 500
)

ggsave(
  filename = file.path("output", "figure", "SupFigS18.pdf"),
  plot = plot_flux_difference("zooplankton") +
    scale_y_continuous(breaks = seq(0, 1, 0.001)),
  height = 8,
  width = 8,
  dpi = 500
)

# Fig 5 ----
message("Generating Fig. 5")
# Preprocessing
tot_null <- timeseries_null |>
  filter(year(sample_week) > 2007) |>
  group_by(station, predator, prey) |>
  summarise(flux = mean(flux, na.rm = TRUE), .groups = "drop")

flux_diff <- station_fluxes |>
  left_join(tot_null, by = c("station", "predator", "prey")) |>
  mutate(
    sig = ifelse(
      (lower - flux > 0),
      "Higher",
      ifelse((upper - flux < 0), "Lower", "Not different")
    ),

    prey = factor(prey, levels = node_data$node_name),
    predator = factor(predator, levels = rev(node_data$node_name)),
    trophic_level = case_when(
      prey %in% node_data$node_name[node_data$trophic_level == 2] ~
        "zooplankton",
      prey %in% node_data$node_name[node_data$trophic_level == 1] ~
        "phytoplankton"
    ),
    rel_mean = mean - flux,
    rel_upper = upper - flux,
    rel_lower = lower - flux,
    rel_null = 0,
    label = ifelse(sig == "Not different", "", round(rel_mean, 2))
  ) |>
  filter(!is.na(flux))

fig5.1 <- ggplot(
  data = filter(flux_diff, trophic_level == "zooplankton"),
  aes(x = prey, y = predator, fill = sig)
) +
  geom_tile(col = "black") +
  coord_fixed() +
  scale_fill_manual(
    values = c(
      "Higher" = "#41E2BA",
      "Lower" = "#2B2D42",
      "Not different" = "white"
    )
  )
fig5.2 <- ggplot(
  data = filter(flux_diff, trophic_level == "phytoplankton"),
  aes(x = prey, y = predator, fill = sig)
) +
  geom_tile(col = "black") +
  coord_fixed() +
  scale_fill_manual(
    values = c(
      "Higher" = "#41E2BA",
      "Lower" = "#2B2D42",
      "Not different" = "white"
    )
  )

#Fig5 <- plot_flux_diff(flux_diff, "zooplankton", ylim = c(-0.0035, 0.0035)) +
#  plot_flux_diff(flux_diff, "phytoplankton", ylim = c(-0.02, 0.02)) +
#  plot_layout(
#    height = c(8, 15),
#    guides = "collect",
#    axes = "collect"
#  )
Fig5 <- (fig5.2 +
  fig5.1 &
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))) +
  plot_layout(guides = "collect", widths = c(12, 8))
ggsave(
  Fig5,
  filename = file.path("output", "figure", "fig5.pdf"),
  width = 12,
  height = 6,
  dpi = 500
)

# Time series analyses -----
# Some timeseries analyses -----
## Temperature ----
temp_data_standard <- temperature |>
  filter(
    station_name == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(year = year(sample_week)) |>
  summarise(temp_avg = mean(temperature, na.rm = T), .groups = "drop") |>
  #  group_by(season) |>
  mutate(
    z = (temp_avg - mean(temp_avg)) / sd(temp_avg),
    col = ifelse(z < 0, "neg", "pos")
  )

# Fig 6 ----
# Some linear regressions
## Fig 6a ----
# Annual Biomass
message("Generating Fig. 6")
annual_biomass <- weekly_biomasses |>
  filter(
    station_name == "BY31 LANDSORTSDJ",
    year(sample_week) > 2007,
    node_name != "Cryptomonadales"
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(year = year(sample_week), season, node_name) |>
  summarise(biomass_avg = mean(biomass, na.rm = T), .groups = "drop") |>
  left_join(node_data, by = join_by(node_name)) |>
  group_by(year, season, trophic_level) |>
  mutate(rel_biomass = biomass_avg / sum(biomass_avg, na.rm = T)) |>
  ungroup() |>
  mutate(
    node_name = factor(node_name, levels = node_data$node_name)
  )
## Fig 6b ----
# Diet overlap
overlap_vs_temp <- inter_overlap_zp |>
  group_by(
    x,
    y,
    sample_week = sample_week.x,
    iso_week = isoweek(sample_week.x),
    year = year(sample_week.x)
  ) |>
  summarise(value = mean(value), .groups = "drop") |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(x, y, year) |>
  summarise(overlap = mean(value), .groups = "drop") |>
  mutate(
    interaction = paste(x, y, sep = "_"),
    interaction = ifelse(
      interaction == "Clupea_Gasterosteus",
      "Gasterosteus_Clupea",
      interaction
    )
  ) |>
  left_join(temp_data_standard, by = join_by(year))

mod_overlap_summary <- run_model_workflow(
  df = overlap_vs_temp,
  group_var = interaction,
  formula_expr = "overlap ~ z",
  filename_suffix = "_overlap_vs_temperature_anomalie.pdf",
  plot_title = "overlap vs temperature anomalie"
)
# primary overlap
primary_vs_temp <- inter_overlap_pp |>
  group_by(
    x,
    y,
    sample_week = sample_week.x,
    iso_week = isoweek(sample_week.x),
    year = year(sample_week.x)
  ) |>
  summarise(value = mean(value), .groups = "drop") |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(x, y, year) |>
  summarise(overlap = mean(value), .groups = "drop") |>
  mutate(
    interaction = paste(x, y, sep = "_"),
    interaction = ifelse(
      interaction == "Clupea_Gasterosteus",
      "Gasterosteus_Clupea",
      interaction
    )
  ) |>
  left_join(temp_data_standard, by = join_by(year))

mod_primary_summary <- run_model_workflow(
  df = primary_vs_temp,
  group_var = interaction,
  formula_expr = "overlap ~ z",
  filename_suffix = "_primary_ overlap_vs_temperature_anomalie.pdf",
  plot_title = "primary overlap vs temperature anomalie"
)
# Based on trophic level (i.e., sum all node from the same trophic level):
tot_biomass <- annual_biomass |>
  group_by(year, type) |>
  summarise(biomass = sum(biomass_avg), .groups = "drop") |>
  left_join(
    temp_data_standard,
    by = join_by(year)
  )
mod_biomass_summary <- run_model_workflow(
  df = tot_biomass,
  group_var = type,
  formula_expr = "biomass ~ z",
  filename_suffix = "_biomass_vs_temperature_anomalie.pdf",
  plot_title = "biomass vs temperature anomalie"
)
## Fig 6c ----
# Predation pressure analyses
#flux_and_temp <- timeseries_fluxes |>
#  left_join(
#    temp_data_standard, # |> filter(season=="Summer"),
#    by = join_by(year)
#  ) |>
#  filter(
#    station == "BY31 LANDSORTSDJ",
#    isoweek(sample_week) %in% 2:51,
#    year(sample_week) %in% 2008:2023
#  ) |>
#  add_season() |>
#  filter(season != "Winter") |>
#  left_join(node_data, by = join_by(predator == node_name)) |>
#  group_by(year = year(sample_week), predator, prey, col, type, z) |>
#  summarise(flux = mean(mean, na.rm = T), .groups = "drop") |>
#  group_by(predator, year, col, type) |>
#  mutate(rel_flux = flux / sum(flux)) |>
#  ungroup()
#tot_consum <-
#  flux_and_temp |>
#  group_by(year, predator, type) |>
#  summarise(flux = mean(flux), .groups = "drop") |>
#  left_join(
#    temp_data_standard, # |> filter(season == "Summer"),
#    by = join_by(year)
#  ) |>
#  mutate(prey = ifelse(type == "fish", "zooplankton", "phytoplankton"))
tot_pred <- pred_pressure_ts |>
  rename("prey" = type, "predation_pressure" = y) |>
  left_join(temp_data_standard, by = join_by(year))

mod_pressure_summary <- run_model_workflow(
  df = tot_pred,
  group_var = prey,
  formula_expr = "predation_pressure ~ z",
  filename_suffix = "_predation_pressure_vs_temperature_anomalie.pdf",
  plot_title = "predation pressure vs temperature anomalie"
)


models_combined <- list(
  predation_pressure = mod_pressure_summary |> rename(group = prey),
  biomass = mod_biomass_summary |> rename(group = type),
  overlap = mod_overlap_summary |> rename(group = interaction),
  primary_overlap = mod_primary_summary |> rename(group = interaction)
) |>
  bind_rows(.id = "model") |>
  ungroup() |>
  filter(term == "slope") |>
  mutate(
    p_adj = p.adjust(p.value, method = "fdr"),
    p_adj_sig = case_when(
      p_adj <= 0.05 ~ "< 0.05",
      p_adj <= 0.1 ~ "< 0.1",
      TRUE ~ "> 0.1"
    ),
    p_adj_sig = factor(p_adj_sig, levels = c("< 0.05", "< 0.1", "> 0.1"))
  ) |>
  select(group, model, term, estimate, r.squared, p.value, p_adj, p_adj_sig)

line_values = c("< 0.05" = 1, "< 0.1" = 2, "> 0.1" = 3)
color_values = c(
  "fish" = "#feedde",
  "phytoplankton" = "#006837",
  "zooplankton" = "#fe9929",
  "Gasterosteus_Clupea" = "#74a9cf",
  "Clupea_Sprattus" = "#ffffd4",
  "Gasterosteus_Sprattus" = "#b30000",
  "ratio" = "black"
)
Fig6.a <- tot_biomass |>
  left_join(
    models_combined |> rename("type" = group) |> filter(model == "biomass"),
    by = join_by(type)
  ) |>
  ggplot(aes(
    x = z,
    y = biomass,
    fill = type,
    col = type,
    linetype = p_adj_sig
  )) +
  geom_point(shape = 21, col = "black", size = 3) +
  stat_smooth(method = "lm", se = F, formula = 'y ~ x') +
  scale_linetype_manual(values = line_values) +
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = color_values) +
  labs(
    x = "Temperature anomalies \n (z-scores)",
    y = "Average biomass from spring to fall \n [g/m2]"
  )
Fig6.b <- tot_pred |>
  left_join(
    models_combined |>
      rename("prey" = group) |>
      filter(model == "predation_pressure"),
    by = join_by(prey)
  ) |>
  ggplot(aes(
    x = z,
    y = predation_pressure,
    fill = prey,
    col = prey,
    linetype = p_adj_sig
  )) +
  geom_point(shape = 21, col = "black", size = 3) +
  stat_smooth(method = "lm", se = F, formula = 'y ~ x') +
  scale_linetype_manual(values = line_values) +
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = color_values) +
  scale_y_continuous(
    sec.axis = sec_axis(
      ~ . / 5,
      name = "Predation pressure \n on zooplankton \n [kJ/g]"
    )
  ) +
  labs(
    x = "Temperature anomalies \n (z-scores)",
    y = "Predation pressure \n on phytoplankton \n [kJ/g]"
  )

Fig6.c <- overlap_vs_temp |>
  left_join(
    models_combined |>
      rename("interaction" = group) |>
      filter(model == "overlap"),
    by = join_by(interaction)
  ) |>
  ggplot(aes(
    x = z,
    y = overlap,
    fill = interaction,
    col = interaction,
    linetype = p_adj_sig
  )) +

  geom_point(shape = 21, col = "black", size = 3) +
  stat_smooth(method = "lm", se = F, formula = 'y ~ x') +
  scale_linetype_manual(values = line_values) +
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = color_values) +
  labs(
    x = "Temperature anomalies \n (z-scores)",
    y = "Diet overlap"
  )
Fig6.d <- primary_vs_temp |>
  left_join(
    models_combined |>
      rename("interaction" = group) |>
      filter(model == "primary_overlap"),
    by = join_by(interaction)
  ) |>
  ggplot(aes(
    x = z,
    y = overlap,
    fill = interaction,
    col = interaction,
    linetype = p_adj_sig
  )) +

  geom_point(shape = 21, col = "black", size = 3) +
  stat_smooth(method = "lm", se = F, formula = 'y ~ x') +
  scale_linetype_manual(values = line_values) +
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = color_values) +
  labs(
    x = "Temperature anomalies \n (z-scores)",
    y = "primary producers overlap"
  )

Fig6 <- Fig6.a /
  Fig6.b /
  Fig6.c /
  Fig6.d +
  plot_layout(guides = "collect", axes = "collect")
ggsave(
  plot = Fig6,
  filename = file.path(
    "output",
    "figure",
    "fig6.pdf"
  ),
  height = 8,
  width = 7
)
# Fig 7 ----
# Correlation analyses
# Biomass anomalies
message("Generating Fig. 7")
biomass_anomalies <- weekly_biomasses |>
  filter(
    station_name == "BY31 LANDSORTSDJ",
    year(sample_week) %in% 2008:2023,
    node_name != "Cryptomonadales",
    node_name != "Cyanobiaceae"
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(node_name, station = station_name, year = year(sample_week)) |>
  summarise(biomass = mean(biomass, na.rm = T), .groups = "drop_last") |>
  mutate(value = (biomass - mean(biomass)) / sd(biomass)) |>
  ungroup() |>
  select(-biomass) |>
  mutate(parameter = "biomass_anomalies")
# Contribution to higher trophic level
contribution <- timeseries_fluxes |>
  filter(
    station == "BY31 LANDSORTSDJ",
    year(sample_week) %in% 2008:2023
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(node_name = prey, station, year = year(sample_week), sample_week) |>
  summarise(flux = sum(mean, na.rm = T), .groups = "drop_last") |>
  summarise(value = mean(flux), .groups = "drop") |>
  left_join(node_data, by = join_by(node_name)) |>
  group_by(type, year) |>
  mutate(value = value / sum(value), parameter = "contribution") |>
  ungroup() |>
  select(node_name, station, year, value, parameter)
# Predation pressure on lower trophic level
predation_pressure_on_lower_trophic_level <- timeseries_fluxes |>
  filter(
    station == "BY31 LANDSORTSDJ",
    year(sample_week) %in% 2008:2023
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(predator, station, year = year(sample_week), sample_week) |>
  summarise(flux = sum(mean, na.rm = T), .groups = "drop_last") |>
  summarise(consumption = mean(flux), .groups = "drop") |>
  mutate(
    prey = ifelse(
      predator %in% c("Clupea", "Sprattus", "Gasterosteus"),
      "zooplankton",
      "phytoplankton"
    )
  ) |>
  left_join(
    weekly_biomasses |>
      filter(
        station_name == "BY31 LANDSORTSDJ",
        year(sample_week) %in% 2008:2023
      ) |>
      add_season() |>
      filter(season != "Winter") |>
      left_join(node_data, join_by(node_name)) |>
      group_by(
        prey = type,
        station_name,
        year = year(sample_week),
        sample_week
      ) |>
      summarise(biomass = sum(biomass, na.rm = T), .groups = "drop_last") |>
      summarise(biomass = mean(biomass), .groups = "drop"),
    by = join_by(year, prey)
  ) |>
  mutate(
    value = consumption / biomass,
    parameter = "predation_pressure_on_lower_trophic_level"
  ) |>
  select(node_name = predator, station, year, value, parameter)
# Predation pressure by higher trophic level
predation_pressure_by_higher_trophic_level <- timeseries_fluxes |>
  filter(
    station == "BY31 LANDSORTSDJ",
    year(sample_week) %in% 2008:2023
  ) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(prey, station, year = year(sample_week), sample_week) |>
  summarise(flux = sum(mean, na.rm = T), .groups = "drop_last") |>
  summarise(contribution = mean(flux), .groups = "drop") |>
  left_join(node_data, by = join_by(prey == node_name)) |>
  mutate(
    predator = ifelse(
      type == "zooplankton",
      "fish",
      "zooplankton"
    )
  ) |>
  left_join(
    weekly_biomasses |>
      filter(
        station_name == "BY31 LANDSORTSDJ",
        year(sample_week) %in% 2008:2023
      ) |>
      add_season() |>
      filter(season != "Winter") |>
      group_by(
        prey = node_name,
        station_name,
        year = year(sample_week),
        sample_week
      ) |>
      summarise(biomass = sum(biomass, na.rm = T), .groups = "drop_last") |>
      summarise(biomass = mean(biomass), .groups = "drop"),
    by = join_by(prey, year)
  ) |>
  mutate(
    value = contribution / biomass,
    parameter = "predation_pressure_by_higher_trophic_level"
  ) |>
  select(node_name = prey, station, year, value, parameter)

# Combine and add temperature anomalies
correlation_df <- bind_rows(
  biomass_anomalies,
  contribution,
  predation_pressure_on_lower_trophic_level,
  predation_pressure_by_higher_trophic_level
) |>
  left_join(temp_data_standard, by = join_by(year))
bg_df <- tibble(
  node_name = 1:(length(node_data$node_name) - 1),
  idx = seq_along(node_name)
) |>
  filter(idx %% 2 == 0)
Fig7 <- correlation_df |>
  group_by(node_name, station, parameter) |>
  summarise(
    cor_test = list(cor.test(value, z, method = "spearman")),
    .groups = "drop"
  ) |>
  mutate(
    result = map(cor_test, tidy)
  ) |>
  unnest(result) |>
  ungroup() |>
  mutate(
    p_adj = p.adjust(p.value, method = "fdr"),
    p_adj_sig = case_when(
      p_adj <= 0.05 ~ "*",
      p_adj <= 0.1 ~ " ",
      TRUE ~ ""
    ),
    p_adj_sig = factor(p_adj_sig, levels = c("*", " ", "")),
    node_name = factor(node_name, levels = node_data$node_name),
    parameter = fct_rev(factor(
      parameter,
      levels = c(
        "biomass_anomalies",
        "contribution",
        "predation_pressure_by_higher_trophic_level",
        "predation_pressure_on_lower_trophic_level"
      )
    ))
  ) |>
  ggplot(aes(
    y = parameter,
    x = node_name,
    fill = estimate,
    label = p_adj_sig
  )) +
  geom_rect(
    data = bg_df,
    aes(
      xmin = node_name - 0.5,
      xmax = node_name + 0.5,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "gray80",
    inherit.aes = FALSE
  ) +
  geom_point(shape = 24, aes(size = abs(estimate))) +
  geom_text() +
  scale_size_continuous(
    range = c(.5, 7),
    limits = c(0, 1),
    breaks = c(0, .4, .8)
  ) +
  scale_fill_gradient2(
    low = "#01665e",
    mid = "white",
    high = "#8c510a",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = seq(-0.8, .8, .4)
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(
  plot = Fig7,
  filename = file.path(
    "output",
    "figure",
    "fig7.pdf"
  ),
  height = 4,
  width = 10
)
message("Everything has run smoothly :)")
