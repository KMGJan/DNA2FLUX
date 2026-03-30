#!/usr/bin/env Rscript

# Wed Dec 17 15:08:39 2025 ------------------------------

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
if (!dir.exists(file.path("output", "figure"))) {
  dir.create(file.path("output", "figure"), recursive = T)
}
if (!dir.exists(file.path("output", "residuals"))) {
  dir.create(file.path("output", "residuals"), recursive = T)
}
if (!dir.exists(file.path("output", "table"))) {
  dir.create(file.path("output", "table"), recursive = T)
}

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
  library(spaa)
  library(broom)
  library(ggrepel)
  library(scales)
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
  # Check if all node have outgoing fluxes and remove the ones that never contribute to higher trophic levels
  group_by(predator, prey) |>
  filter(
    sum(flux, na.rm = T) > 0,
    # only keep year 2008-2023 and weeks 2-51
    year(sample_week) > 2007,
    isoweek(sample_week) %in% 2:51
  ) |>
  ungroup()

timeseries_fluxes <- readAndArrange(
  file = file.path("data", "analyses", "timeseries_fluxes.csv")
) |>
  # Check if all node have outgoing fluxes and remove the ones that never contribute to higher trophic levels
  group_by(predator, prey) |>
  filter(
    sum(mean, na.rm = T) > 0,
    # only keep year 2008-2023 and weeks 2-51
    year(sample_week) > 2007,
    isoweek(sample_week) %in% 2:51
  ) |>
  ungroup()

isoweek_fluxes <- readAndArrange(
  file = file.path("data", "analyses", "isoweek_fluxes.csv"),
  arrange = FALSE
) |>
  filter(iso_week %in% 2:51)

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
  # Check if all node are present and remove the ones that never occur
  filter(
    # Ensure that only the
    node_name %in%
      unique(timeseries_fluxes$prey) |
      node_name %in% unique(timeseries_fluxes$predator),
    # only keep BY31
    station_name == unique(timeseries_fluxes$station),
    # only keep year 2008-2023 and weeks 2-51
    year(sample_week) > 2007,
    isoweek(sample_week) %in% 2:51
  )

temperature <- readAndArrange(
  file = file.path("data", "processed", "interpolation", "temperature.csv")
) |>
  filter(
    # only keep year 2008-2023 and weeks 2-51
    year(sample_week) %in% 2008:2023,
    isoweek(sample_week) %in% 2:51,
    # only keep BY31
    station_name == "BY31 LANDSORTSDJ"
  )
# Set the color scheme for plotting
color_mapping <- setNames(node_data$color, node_data$node_name)

# Start the analyses ---------
# Fig. S6: input data ----
message("Fig. S6 Input data")
# Create the seasonal input data
seasonal_input <- weekly_biomasses |>
  left_join(node_data, by = join_by(node_name)) |>
  filter(type != "fish") |>
  mutate(
    node_name = factor(node_name, levels = node_data$node_name),
    type = fct_rev(factor(
      type,
      levels = c("fish", "zooplankton", "phytoplankton")
    ))
  ) |>
  # Group by iso week and get the average +- sd of the biomass for each node
  group_by(type, node_name, iso_week = isoweek(sample_week)) |>
  summarise(
    y = mean(biomass, na.rm = T),
    ymin = pmax(0, y - sd(biomass, na.rm = T)),
    ymax = y + sd(biomass, na.rm = T),
    .groups = "drop"
  )
# Add empty rows for fish
fish_placeholder <- tibble(
  type = factor("fish", levels = c("fish", "zooplankton", "phytoplankton")),
  node_name = "Clupea",
  iso_week = 2,
  y = NA_real_,
  ymin = 0,
  ymax = 0
)
# Fish relative biomass over years
FigS6h <-
  weekly_biomasses |>
  left_join(node_data, by = join_by(node_name)) |>
  filter(type == "fish") |>
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
# Plankton biomass over seasons
FigS6fg <-
  seasonal_input |>
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
# Temperature over seasons
FigS6e <-
  temperature |>
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
# Combine plots and save
FigS6efgh <- FigS6e /
  FigS6fg /
  FigS6h +
  plot_layout(heights = c(1, 2, 1), guides = "collect", axes = "collect") &
  guides(colour = guide_legend(ncol = 1))
ggsave(
  plot = FigS6efgh,
  filename = file.path("output", "figure", "SupFigS6efgh.pdf"),
  width = 5,
  height = 8
)

# Interannual dynamics
# Plot background dataset
bg_df <- tibble(
  year = 2008:2023,
  idx = seq_along(year)
) |>
  filter(idx %% 2 == 0)

# Aggreagated trophic level average biomass over years
FigS6bcd <-
  weekly_biomasses |>
  left_join(node_data, by = join_by(node_name)) |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(type, year = year(sample_week), sample_week) |>
  summarise(biomass = sum(biomass), .groups = "drop_last") |>
  summarise(biomass = mean(biomass, na.rm = TRUE), .groups = "drop_last") |>
  # compute interannual average for each trophic level
  mutate(
    avg = mean(biomass),
    type = factor(type, levels = c("phytoplankton", "zooplankton", "fish"))
  ) |>
  ungroup() |>
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
  # Lollipop segments
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
  labs(x = NULL, y = "Biomass")

# Aggreagated temperature over years
FigS6a <-
  temperature |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(year = year(sample_week), type = "temperature") |>
  summarise(temp_avg = mean(temperature, na.rm = T), .groups = "drop") |>
  # Compute average temperature overall
  mutate(avg = mean(temp_avg)) |>
  # Plot
  ggplot(aes(x = year, y = temp_avg)) +
  # Background rectangles
  geom_rect(
    data = bg_df,
    aes(xmin = year - 0.5, xmax = year + 0.5, ymin = -Inf, ymax = Inf),
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
FigS6abcd <- FigS6a /
  FigS6bcd +
  plot_layout(heights = c(1, 3), guides = "collect", axes = "collect")
ggsave(
  plot = FigS6abcd,
  filename = file.path("output", "figure", "SupFigS1abcd.pdf"),
  width = 8,
  height = 8
)

message("Fig. S6 plotted")
# Add the stats ....
df_test_temperature <-
  temperature |>
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
# Fig. 3: PCoA --------
## Produce the dataframes ----
### Zooplankton --> fish relative flux  ----
rel_flux_zp <-
  timeseries_fluxes |>
  filter() |>
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

### Primary-producer --> fish relative flux ----
rel_flux_pp <- timeseries_fluxes |>
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
    rel_flux_zp |>
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

### Ambient zooplankton and phytoplankton relative contribution ----
rel_biomass <-
  weekly_biomasses |>
  mutate(
    TL = case_when(
      node_name %in% node_data$node_name[node_data$trophic_level == 3] ~ "fish",
      node_name %in% node_data$node_name[node_data$trophic_level == 2] ~
        "zooplankton",
      node_name %in% node_data$node_name[node_data$trophic_level == 1] ~
        "phytoplankton"
    )
  ) |>
  group_by(sample_week, TL) |>
  filter(sum(biomass) > 0) |>
  mutate(rel_biomass = biomass / sum(biomass)) |>
  ungroup()

## From datafram to matrices ----
mat_zp <- rel_flux_zp |>
  filter(type == "fish") |>
  mutate(type = "selectivity") |>
  pivot_wider(
    names_from = prey,
    values_from = rel_flux,
    values_fill = 0
  ) |>
  na.omit()
mat_pp <- rel_flux_pp |>
  mutate(type = "selectivity") |>
  pivot_wider(
    names_from = primary_producer,
    values_from = rel_primary_production,
    values_fill = 0
  ) |>
  na.omit()
### Add the ambient primary producers and zooplankton community
ambient_zp <-
  rel_biomass |>
  filter(TL == "zooplankton") |>
  mutate(predator = "ambient", type = "ambient") |>
  select(
    node_name,
    sample_week,
    predator,
    station = "station_name",
    type,
    rel_biomass
  ) |>
  pivot_wider(
    names_from = node_name,
    values_from = rel_biomass,
    values_fill = 0
  )
ambient_pp <-
  rel_biomass |>
  filter(TL == "phytoplankton") |>
  mutate(predator = "ambient", type = "ambient") |>
  select(
    node_name,
    sample_week,
    predator,
    station = "station_name",
    type,
    rel_biomass
  ) |>
  pivot_wider(
    names_from = node_name,
    values_from = rel_biomass,
    values_fill = 0
  )

### combine these to the respective matrices
mat_zp_complete <-
  mat_zp |>
  bind_rows(ambient_zp)
mat_pp_complete <-
  mat_pp |>
  bind_rows(ambient_pp)
# and one matrix for phytoplankton --> zooplankton
mat_pp2zp <- rel_flux_zp |>
  filter(
    type == "zooplankton" #,
    #For faster computing time, only keep even isoweek to reduce the size of the dataset
    #isoweek(sample_week) %% 2 == 0
  ) |>
  mutate(type = "selectivity") |>
  pivot_wider(
    names_from = prey,
    values_from = rel_flux,
    values_fill = 0
  ) |>
  na.omit() |>
  bind_rows(ambient_pp)
### Run the PCoA -----
matrices <- list(
  pp = mat_pp_complete,
  zp = mat_zp_complete,
  pp2zp = mat_pp2zp
)
# Because it is long, make it only run if it hasn't already...
if (!file.exists(file.path("output", "table", "pcoa.rds"))) {
  message("The 3 PCoAs are running in parrallel... this takes about 10 min...")
  suppressPackageStartupMessages(library(furrr))
  plan(multisession)
  pcoa_workflow <- future_map(
    matrices,
    \(mat) run_pcoa(mat),
    .options = furrr_options(seed = TRUE)
  )
  pcoa_workflow |> saveRDS(file.path("output", "table", "pcoa.rds"))
} else {
  pcoa_workflow <- readRDS(file.path("output", "table", "pcoa.rds"))
}
## primary producers --> fish
pcoa_df_pp <- pcoa_workflow$pp$pcoa_df
eig_pp <- pcoa_workflow$pp$eig
## zooplankton --> fish
pcoa_df_zp <- pcoa_workflow$zp$pcoa_df
eig_zp <- pcoa_workflow$zp$eig
# primary producers --> zooplankton
pcoa_df_pp2zp <- pcoa_workflow$pp2zp$pcoa_df
eig_pp2zp <- pcoa_workflow$pp2zp$eig

## Visualisation ----
# bind envfit with type
envfit <- bind_rows(
  pcoa_workflow$zp$envfit |> mutate(source = "zooplankton"),
  pcoa_workflow$pp$envfit |> mutate(source = "phytoplankton")
)
# Plot
message("Fig. 3 PCoA")

# Mon Mar 30 10:36:01 2026 ------------------------------

Fig3 <- bind_rows(
  pcoa_df_zp |> mutate(source = "zooplankton"),
  pcoa_df_pp |> mutate(source = "phytoplankton")
) |>
  mutate(
    source = factor(source, levels = c("zooplankton", "phytoplankton")),
    predator = factor(predator, levels = c('ambient', node_data$node_name))
  ) |>
  group_by(iso_week, predator, season, type, source) |>
  summarise(
    avg1 = mean(Axis1, na.rm = TRUE),
    avg2 = mean(Axis2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(
    x = avg1,
    y = avg2,
    fill = iso_week,
    col = iso_week
  )) +
  geom_segment(
    data = envfit |>
      mutate(
        source = factor(source, levels = c("zooplankton", "phytoplankton"))
      ),
    aes(x = Axis1, y = Axis2, xend = 0, yend = 0),
    arrow = arrow(length = unit(2, 'mm'), ends = "first"),
    inherit.aes = FALSE
  ) +
  geom_text_repel(
    data = envfit |>
      mutate(
        source = factor(source, levels = c("zooplankton", "phytoplankton"))
      ),
    aes(x = Axis1, y = Axis2, label = prey),
    inherit.aes = FALSE,
    size = 3,
    seed = 100
  ) +
  geom_point(size = 2, col = "black", shape = 21) +
  facet_grid(source ~ predator) +
  coord_fixed() +
  scale_fill_gradientn(
    colors = c("#e66101", "white", "#5e3c99"),
    values = scales::rescale(c(2, 22, 48)),
    limits = c(1, 52),
    breaks = c(12, 24, 36, 48),
    labels = c("Mar", "Jun", "Sep", "Dec")
  ) +
  labs(x = "PCOA1", y = "PCOA2")
# Save
ggsave(
  plot = Fig3,
  filename = file.path("output", "figure", "fig3.pdf"),
  dpi = 500,
  width = 7.5,
  height = 4
)
message(
  "PCoA1 explains ",
  round(eig_zp[c(1)] * 100, 1),
  "% and PCoA2 ",
  round(eig_zp[c(2)] * 100, 1),
  "% of the variation in the diet composition ordination"
)
message(
  "PCoA1 explains ",
  round(eig_pp[c(1)] * 100, 1),
  "% and PCoA2 ",
  round(eig_pp[c(2)] * 100, 1),
  "% of the variation in the primary producer sources ordination"
)

rm(eig_pp, eig_zp, envfit)
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
    formula = pcoa_workflow$zp$bray ~ predator + factor(iso_week),
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
    formula = pcoa_workflow$pp$bray ~ predator + factor(iso_week),
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
# Fig S5: Diet overlap ----
message("Fig. S5, diet overlap over time")
inter_overlap_zp <-
  mat_zp |>
  filter(isoweek(sample_week) %in% 2:51) |>
  mutate(
    group_id = paste(year(sample_week), isoweek(sample_week), sep = "_")
  ) |>
  group_by(group_id) |>
  group_split() |>
  map_dfr(
    ~ interspecific_overlap(.x, 5:12)
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
    ~ interspecific_overlap(.x, 5:16)
  ) |>
  mutate(facet = "interspecific") |>
  filter(x != y)

FigS5bd <-
  inter_overlap_zp |>
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
FigS5ac <-
  inter_overlap_zp |>
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
FigS5 <- FigS5ac + FigS5bd + plot_layout(guides = "collect")
ggsave(
  plot = FigS5,
  filename = file.path("output", "figure", "SupFigS5.pdf"),
  width = 9,
  height = 4
)
rm(FigS5ac, FigS5bd)

# Fig. 2: Food web topology -----
# Save the PCoA1 loadings for plotting the foodweb with and without selectivity:
message("Fig. 2 Food web topology")
position_zp <- make_seasonal_position(pcoa_df_pp2zp)
position_fish <- make_seasonal_position(pcoa_df_zp) |>
  mutate(position_x = 1 - position_x)
position_pp <-
  pcoa_workflow$pp2zp$envfit |>
  mutate(predator = prey, position_x = scales::rescale(Axis1, c(0, 1))) |>
  select(predator, position_x) |>
  mutate(
    position_y = case_when(
      predator == "Peridiniales" ~ 1.02,
      predator == "Thalassiosirales" ~ .98,
      .default = 1
    )
  ) |>
  cross_join(tibble(season = c("Spring", "Summer", "Fall", "Winter"))) |>
  cross_join(tibble(type = c("ambient", "selectivity")))

# Mon Mar 30 11:41:19 2026 ------------------------------
# Add some jitter for the fish and zooplankton in the neutral model
set.seed(10)
position_data <-
  bind_rows(
    position_pp,
    position_zp |>
      filter(type == 'ambient') |>
      select(-predator) |>
      cross_join(tibble(
        predator = node_data$node_name[node_data$type == 'zooplankton']
      )) |>
      group_by(type, season) |>
      mutate(
        position_x = jitter(position_x, amount = .05),
        position_y = 2,
        position_y = jitter(position_y, amount = .15)
      ) |>
      ungroup() |>
      bind_rows(
        position_zp |>
          filter(type != 'ambient') |>
          mutate(position_y = 2)
      ),
    position_fish |>
      filter(type == 'ambient') |>
      select(-predator) |>
      cross_join(tibble(
        predator = node_data$node_name[node_data$type == 'fish']
      )) |>
      group_by(type, season) |>
      mutate(
        position_x = jitter(position_x, amount = .05),
        position_y = 3,
        position_y = jitter(position_y, amount = .15)
      ) |>
      ungroup() |>
      bind_rows(
        position_fish |>
          filter(type != 'ambient') |>
          mutate(position_y = 3)
      )
  ) |>
  mutate(
    season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter")),
    predator = factor(predator, levels = node_data$node_name)
  )

# Run the function food_web_fig2 for each season and selecitivity combination
seasons <- c("Spring", "Summer", "Fall", "Winter")
plots <-
  c(
    lapply(seasons, \(x) food_web_fig2(x, selectivity = FALSE)),
    lapply(seasons, \(x) food_web_fig2(x, selectivity = TRUE))
  )
Fig2 <- wrap_plots(plots, guides = "collect", nrow = 2)
ggsave(
  plot = Fig2,
  filename = file.path("output", "figure", "Fig2.pdf"),
  width = 13,
  height = 10
)
# FigS4: Relative contribution of prey for each trophic level, per season ----
## With selectivity
rel_contrib_select <- timeseries_fluxes |>
  filter(
    station == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  add_season() |>
  group_by(predator, prey, station, season, year) |>
  summarise(flux = mean(mean, na.rm = TRUE), .groups = "drop_last") |>
  summarise(flux_avg = mean(flux, na.rm = TRUE), .groups = "drop") |>
  mutate(
    trophic_level = case_when(
      predator %in% node_data$node_name[node_data$trophic_level == 1] ~
        "phytoplankton",
      predator %in% node_data$node_name[node_data$trophic_level == 2] ~
        "zooplankton",
      .default = "fish"
    )
  ) |>
  group_by(trophic_level, season, prey) |>
  summarise(flux = sum(flux_avg), .groups = "drop_last") |>
  mutate(
    relFlux = 100 * flux / sum(flux)
  ) |>
  ungroup() |>
  rbind(
    timeseries_fluxes |>
      filter(
        station == "BY31 LANDSORTSDJ",
        isoweek(sample_week) %in% 2:51,
        year(sample_week) %in% 2008:2023
      ) |>

      group_by(predator, prey, station, year) |>
      summarise(flux = mean(mean, na.rm = TRUE), .groups = "drop_last") |>
      summarise(flux_avg = mean(flux, na.rm = TRUE), .groups = "drop") |>
      mutate(
        trophic_level = case_when(
          predator %in% node_data$node_name[node_data$trophic_level == 1] ~
            "phytoplankton",
          predator %in% node_data$node_name[node_data$trophic_level == 2] ~
            "zooplankton",
          .default = "fish"
        )
      ) |>
      group_by(trophic_level, prey) |>
      summarise(flux = sum(flux_avg), .groups = "drop_last") |>
      mutate(
        relFlux = 100 * flux / sum(flux),
        prey = factor(prey, levels = node_data$node_name),
        season = "annual"
      )
  ) |>
  ungroup()
## Withoutselectivity
rel_contrib_Noselect <-
  timeseries_null |>
  filter(
    station == "BY31 LANDSORTSDJ",
    isoweek(sample_week) %in% 2:51,
    year(sample_week) %in% 2008:2023
  ) |>
  add_season() |>
  group_by(predator, prey, station, season, year = year(sample_week)) |>
  summarise(flux = mean(flux, na.rm = TRUE), .groups = "drop_last") |>
  summarise(flux_avg = mean(flux, na.rm = TRUE), .groups = "drop") |>
  mutate(
    trophic_level = case_when(
      predator %in% node_data$node_name[node_data$trophic_level == 1] ~
        "phytoplankton",
      predator %in% node_data$node_name[node_data$trophic_level == 2] ~
        "zooplankton",
      .default = "fish"
    )
  ) |>
  group_by(trophic_level, season, prey) |>
  summarise(flux = sum(flux_avg), .groups = "drop_last") |>
  mutate(
    relFlux = 100 * flux / sum(flux)
  ) |>
  ungroup() |>
  rbind(
    timeseries_null |>
      filter(
        station == "BY31 LANDSORTSDJ",
        isoweek(sample_week) %in% 2:51,
        year(sample_week) %in% 2008:2023
      ) |>
      group_by(predator, prey, station, year = year(sample_week)) |>
      summarise(flux = mean(flux, na.rm = TRUE), .groups = "drop_last") |>
      summarise(flux_avg = mean(flux, na.rm = TRUE), .groups = "drop") |>
      mutate(
        trophic_level = case_when(
          predator %in% node_data$node_name[node_data$trophic_level == 1] ~
            "phytoplankton",
          predator %in% node_data$node_name[node_data$trophic_level == 2] ~
            "zooplankton",
          .default = "fish"
        )
      ) |>
      group_by(trophic_level, prey) |>
      summarise(flux = sum(flux_avg), .groups = "drop_last") |>
      mutate(
        relFlux = 100 * flux / sum(flux),
        season = "annual"
      ) |>
      ungroup()
  )
# Plot
FigS4 <- rel_contrib_Noselect |>
  mutate(selectivity = F) |>
  bind_rows(rel_contrib_select |> mutate(selectivity = T)) |>
  mutate(
    trophic_level = factor(trophic_level, levels = c("zooplankton", "fish")),
    season = factor(
      season,
      levels = c("annual", "Spring", "Summer", "Fall", "Winter")
    ),
    prey = factor(prey, levels = node_data$node_name)
  ) |>
  ggplot(aes(x = season, y = relFlux, fill = prey)) +
  geom_bar(stat = "identity") +
  facet_grid(selectivity ~ trophic_level) +
  scale_fill_manual(values = color_mapping) +
  theme(legend.position = "bottom")
ggsave(
  plot = FigS4,
  filename = file.path("output", "figure", "SupFigS4.pdf"),
  width = 8,
  height = 6
)

# Fig. 1: first and last food web -----
message("Fig. 1 First and last food web")
position_zp <- make_position_fig1(pcoa_df_pp2zp)
position_fish <- make_position_fig1(pcoa_df_zp) |>
  mutate(position_x = 1 - position_x)
position_pp <-
  pcoa_workflow$pp2zp$envfit |>
  mutate(predator = prey, position_x = scales::rescale(Axis1, c(0, 1))) |>
  select(predator, position_x) |>
  cross_join(tibble(sample_week = unique(position_zp$sample_week))) |>
  cross_join(tibble(type = c("ambient", "selectivity")))
# Mon Mar 30 13:38:04 2026 ------------------------------
# Adding jitter for zooplankton and fish under neutral food web
set.seed(10)
position_data <- bind_rows(
  position_pp |> mutate(position_y = 1),
  position_zp |>
    filter(type == 'ambient') |>
    select(-predator) |>
    cross_join(tibble(
      predator = node_data$node_name[node_data$type == 'zooplankton']
    )) |>
    mutate(
      position_x = jitter(position_x, amount = .05),
      position_y = 2,
      position_y = jitter(position_y, amount = .1)
    ) |>
    bind_rows(
      position_zp |>
        filter(type != 'ambient') |>
        mutate(position_y = 2)
    ),
  position_fish |>
    filter(type == 'ambient') |>
    select(-predator) |>
    cross_join(tibble(
      predator = node_data$node_name[node_data$type == 'fish']
    )) |>
    mutate(
      position_x = jitter(position_x, amount = .05),
      position_y = 3,
      position_y = jitter(position_y, amount = .1)
    ) |>
    bind_rows(
      position_fish |>
        filter(type != 'ambient') |>
        mutate(position_y = 3)
    )
)

plots <-
  c(
    lapply(c(T, F), \(x) food_web_fig1(x, selectivity = FALSE)),
    lapply(c(T, F), \(x) food_web_fig1(x, selectivity = TRUE))
  )
fig1 <- wrap_plots(plots, guides = "collect", nrow = 2)
ggsave(
  plot = fig1,
  filename = file.path("output", "figure", "Fig1.pdf"),
  width = 5,
  height = 5
)
# Fig. 4: Network Analysis -----
message("Fig. 4 some network analyses")
network_metrics_df <-
  timeseries_fluxes |>
  group_by(sample_week) |>
  group_map(
    ~ {
      mat <- .x |>
        select(predator, prey, mean) |>
        pivot_wider(values_from = mean, names_from = prey, values_fill = 0) |>
        column_to_rownames("predator") |>
        as.matrix() |>
        t()
      tibble(
        sample_week = as.Date(.y$sample_week),
        connectance = lw(mat, parameter = "connectance"),
        vulnerability = lw(mat, parameter = "vulnerability"),
        generality = lw(mat, parameter = "generality")
      )
    }
  ) |>
  bind_rows() |>
  add_season()
network_metrics_null <-
  timeseries_null |>
  group_by(sample_week) |>
  group_map(
    ~ {
      mat <- .x |>
        select(predator, prey, flux) |>
        pivot_wider(values_from = flux, names_from = prey, values_fill = 0) |>
        column_to_rownames("predator") |>
        as.matrix() |>
        t()
      tibble(
        sample_week = as.Date(.y$sample_week),
        connectance = lw(mat, parameter = "connectance"),
        vulnerability = lw(mat, parameter = "vulnerability"),
        generality = lw(mat, parameter = "generality")
      )
    }
  ) |>
  bind_rows() |>
  add_season()

Fig4 <-
  network_metrics_df |>
  pivot_longer(2:4, values_to = "values", names_to = "parameter") |>
  mutate(selectivity = "selectivity") |>
  rbind(
    network_metrics_null |>
      pivot_longer(2:4, values_to = "values", names_to = "parameter") |>
      mutate(selectivity = "noSelectivity")
  ) |>
  group_by(iso_week, parameter, selectivity) |>
  summarise(avg = mean(values), SD = sd(values), .groups = "drop") |>
  mutate(parameter = paste("link weighted", parameter)) |>
  ggplot(aes(
    x = iso_week,
    y = avg,
    ymin = avg - SD,
    ymax = avg + SD,
    col = selectivity,
    fill = selectivity
  )) +
  geom_ribbon(alpha = .5, col = NA) +
  geom_line(linewidth = 1) +
  facet_grid(parameter ~ ., scales = "free", switch = "y") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.text = element_text(color = "black")
  ) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = c("#466995", "#9F8B7C")) +
  scale_color_manual(values = c("#466995", "#9F8B7C"))
ggsave(
  plot = Fig4,
  filename = file.path("output", "figure", "fig4.pdf"),
  dpi = 500,
  width = 5,
  height = 7
)

# Fig. 5: Warming impacts on diet overlap ----
message("Fig. 5 how warming impacts diet overlap?")
overlap_vs_temp <-
  inter_overlap_zp |>
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
  left_join(df_test_temperature, by = join_by(year))

mod_overlap_summary <- run_model_workflow(
  df = overlap_vs_temp,
  group_var = interaction,
  formula_expr = "overlap ~ z",
  filename_suffix = "_overlap_vs_temperature_anomalie.pdf",
  plot_title = "overlap vs temperature anomalie"
) |>
  filter(term == "slope") |>
  ungroup() |>
  mutate(
    p_adj = p.adjust(p.value, method = "fdr"),
    p_adj_sig = case_when(
      p_adj <= 0.05 ~ "< 0.05",
      p_adj <= 0.1 ~ "< 0.1",
      TRUE ~ "> 0.1"
    ),
    p_adj_sig = factor(p_adj_sig, levels = c("< 0.05", "< 0.1", "> 0.1"))
  ) |>
  select(
    interaction,
    model,
    term,
    estimate,
    r.squared,
    p.value,
    p_adj,
    p_adj_sig
  )

Fig5a <- overlap_vs_temp |>
  left_join(
    mod_overlap_summary,
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
  scale_linetype_manual(values = c("< 0.05" = 1, "< 0.1" = 2, "> 0.1" = NA)) +
  scale_color_manual(
    values = c(
      "Gasterosteus_Clupea" = "#74a9cf",
      "Clupea_Sprattus" = "#ffffd4",
      "Gasterosteus_Sprattus" = "#b30000"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Gasterosteus_Clupea" = "#74a9cf",
      "Clupea_Sprattus" = "#ffffd4",
      "Gasterosteus_Sprattus" = "#b30000"
    )
  ) +
  labs(
    x = "Temperature anomalies \n (z-scores)",
    y = "Diet overlap"
  )

# Biomass anomalies
biomass_anomalies <-
  weekly_biomasses |>
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
contribution <-
  timeseries_fluxes |>
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
  group_by(node_name, station, parameter = "contribution_anomalies") |>
  mutate(value = (value - mean(value)) / sd(value)) |>
  ungroup() |>
  select(node_name, station, year, value, parameter)
# Combine and add temperature anomalies
correlation_results <-
  bind_rows(
    biomass_anomalies,
    contribution
  ) |>
  left_join(df_test_temperature, by = join_by(year)) |>
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
      p_adj <= 0.1 ~ ".",
      TRUE ~ ""
    ),
    p_adj_sig = factor(p_adj_sig, levels = c("*", ".", "")),
    node_name = factor(node_name, levels = node_data$node_name),
    parameter = fct_rev(factor(
      parameter,
      levels = c(
        "biomass_anomalies",
        "contribution_anomalies"
      )
    ))
  )

Fig5b <-
  bind_rows(
    biomass_anomalies,
    contribution
  ) |>
  left_join(df_test_temperature, by = join_by(year)) |>
  left_join(correlation_results, by = join_by(node_name, station, parameter)) |>
  mutate(
    TL = case_when(
      node_name %in% node_data$node_name[node_data$trophic_level == 3] ~ "fish",
      node_name %in% node_data$node_name[node_data$trophic_level == 2] ~
        "zooplankton",
      node_name %in% node_data$node_name[node_data$trophic_level == 1] ~
        "phytoplankton"
    ),
    TL = factor(TL, levels = c("phytoplankton", "zooplankton", "fish"))
  ) |>
  filter(p_adj <= .05) |>
  ggplot(aes(
    x = z,
    y = value,
    col = node_name,
    fill = node_name
  )) +
  geom_point(col = 1, shape = 21) +
  geom_smooth(method = "lm", se = F, formula = 'y ~ x') +
  facet_grid(TL ~ parameter) +
  scale_fill_manual(values = color_mapping) +
  scale_color_manual(values = color_mapping)
Fig5 <- Fig5a / Fig5b
ggsave(
  plot = Fig5,
  filename = file.path("output", "figure", "fig5.pdf"),
  dpi = 500,
  width = 7,
  height = 8
)
# Fig S7 - S8: all correlations -----
FigS7 <-
  biomass_anomalies |>
  left_join(df_test_temperature, by = join_by(year)) |>
  left_join(correlation_results, by = join_by(node_name, station, parameter)) |>
  mutate(
    TL = case_when(
      node_name %in% node_data$node_name[node_data$trophic_level == 3] ~ "fish",
      node_name %in% node_data$node_name[node_data$trophic_level == 2] ~
        "zooplankton",
      node_name %in% node_data$node_name[node_data$trophic_level == 1] ~
        "phytoplankton"
    ),
    TL = factor(TL, levels = c("phytoplankton", "zooplankton", "fish")),
    node_name = factor(node_name, levels = node_data$node_name)
  ) |>
  ggplot(aes(
    x = z,
    y = value,
    col = node_name,
    fill = node_name
  )) +
  geom_point(col = 1, shape = 21) +
  #geom_smooth(method = "lm", se = F, formula = 'y ~ x') +
  facet_wrap(. ~ node_name) +
  scale_fill_manual(values = color_mapping) +
  #scale_color_manual(values = color_mapping) +
  geom_text(
    data = biomass_anomalies |>
      left_join(df_test_temperature, by = join_by(year)) |>
      left_join(
        correlation_results,
        by = join_by(node_name, station, parameter)
      ) |>
      mutate(
        TL = case_when(
          node_name %in% node_data$node_name[node_data$trophic_level == 3] ~
            "fish",
          node_name %in% node_data$node_name[node_data$trophic_level == 2] ~
            "zooplankton",
          node_name %in% node_data$node_name[node_data$trophic_level == 1] ~
            "phytoplankton"
        ),
        TL = factor(TL, levels = c("phytoplankton", "zooplankton", "fish")),
        node_name = factor(node_name, levels = node_data$node_name)
      ) |>
      select(node_name, parameter, estimate, p_adj) |>
      unique(),
    mapping = aes(
      x = 0,
      y = 2.2,
      label = paste0("rho = ", round(estimate, 2), "; P = ", round(p_adj, 3))
    ),
    inherit.aes = FALSE,
    size = 2.7
  ) +
  theme(legend.position = "none") +
  labs(
    x = "Temperature anomalies (z-scores)",
    y = "Biomass anomalies (z-scores)"
  )
ggsave(
  plot = FigS7,
  filename = file.path("output", "figure", "SupFigS7.pdf"),
  dpi = 500,
  width = 7,
  height = 8
)
FigS8 <-
  contribution |>
  left_join(df_test_temperature, by = join_by(year)) |>
  left_join(correlation_results, by = join_by(node_name, station, parameter)) |>
  mutate(
    TL = case_when(
      node_name %in% node_data$node_name[node_data$trophic_level == 3] ~ "fish",
      node_name %in% node_data$node_name[node_data$trophic_level == 2] ~
        "zooplankton",
      node_name %in% node_data$node_name[node_data$trophic_level == 1] ~
        "phytoplankton"
    ),
    TL = factor(TL, levels = c("phytoplankton", "zooplankton", "fish")),
    node_name = factor(node_name, levels = node_data$node_name)
  ) |>
  ggplot(aes(
    x = z,
    y = value,
    col = node_name,
    fill = node_name
  )) +
  geom_point(col = 1, shape = 21) +
  #geom_smooth(method = "lm", se = F, formula = 'y ~ x') +
  facet_wrap(. ~ node_name) +
  scale_fill_manual(values = color_mapping) +
  #scale_color_manual(values = color_mapping) +
  geom_text(
    data = contribution |>
      left_join(df_test_temperature, by = join_by(year)) |>
      left_join(
        correlation_results,
        by = join_by(node_name, station, parameter)
      ) |>
      mutate(
        TL = case_when(
          node_name %in% node_data$node_name[node_data$trophic_level == 3] ~
            "fish",
          node_name %in% node_data$node_name[node_data$trophic_level == 2] ~
            "zooplankton",
          node_name %in% node_data$node_name[node_data$trophic_level == 1] ~
            "phytoplankton"
        ),
        TL = factor(TL, levels = c("phytoplankton", "zooplankton", "fish")),
        node_name = factor(node_name, levels = node_data$node_name)
      ) |>
      select(node_name, parameter, estimate, p_adj) |>
      unique(),
    mapping = aes(
      x = 0,
      y = 2.2,
      label = paste0("rho = ", round(estimate, 2), "; P = ", round(p_adj, 3))
    ),
    inherit.aes = FALSE,
    size = 2.7
  ) +
  theme(legend.position = "none") +
  labs(
    x = "Temperature anomalies (z-scores)",
    y = "Relative outgoing flux anomalies (z-scores)"
  )
ggsave(
  plot = FigS8,
  filename = file.path("output", "figure", "SupFigS8.pdf"),
  dpi = 500,
  width = 7,
  height = 7
)
# Fig S1: Biomass, outgoing fluxes and predation pressure dynamics ----
message("Fig. S1 Biomass, outgoing fluxes and predation pressure dynamics")
figs1_df <- timeseries_fluxes |>
  left_join(node_data, by = c(prey = "node_name")) |>
  left_join(weekly_biomasses, by = join_by(sample_week, prey == node_name)) |>
  group_by(
    trophic_level = case_when(
      trophic_level == 1 ~ "phytoplankton",
      trophic_level == 2 ~ "zooplankton",
      trophic_level == 3 ~ "fish"
    ),
    iso_week,
    sample_week
  ) |>
  summarise(
    outgoing_flux = sum(mean),
    standing_biomass = sum(biomass),
    predation_pressure = outgoing_flux /
      standing_biomass,
    .groups = "drop"
  )

FigS1a <- figs1_df |>
  pivot_longer(4:6, names_to = "parameter") |>
  mutate(
    value = ifelse(
      parameter == "outgoing_flux" & trophic_level == "zooplankton",
      value * 10,
      value
    ),
    parameter = factor(
      parameter,
      levels = c("standing_biomass", "outgoing_flux", "predation_pressure")
    )
  ) |>
  group_by(trophic_level, iso_week, parameter) |>
  summarise(avg = mean(value), sd = sd(value), .groups = "drop") |>
  ggplot(aes(
    x = iso_week,
    y = avg,
    ymin = avg - sd,
    ymax = avg + sd,
    fill = trophic_level,
    col = trophic_level
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
FigS1b <- figs1_df |>
  pivot_longer(4:6, names_to = "parameter") |>
  mutate(
    value = ifelse(
      parameter == "outgoing_flux" & trophic_level == "zooplankton",
      value * 10,
      value
    ),
    parameter = factor(
      parameter,
      levels = c("standing_biomass", "outgoing_flux", "predation_pressure")
    )
  ) |>
  filter(iso_week %in% 11:48) |>
  group_by(trophic_level, year = year(sample_week), parameter) |>
  summarise(avg = mean(value), .groups = "drop") |>
  ggplot(aes(x = year, y = avg, fill = trophic_level, col = trophic_level)) +
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
FigS1 <- FigS1a + FigS1b
ggsave(
  filename = file.path("output", "figure", "SupFigS1.pdf"),
  plot = FigS1,
  height = 7,
  width = 6,
  dpi = 500
)
tot_flux <- timeseries_fluxes |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(
    trophic_level = ifelse(
      predator %in% node_data$node_name[node_data$trophic_level == 2],
      "zooplankton",
      "fish"
    ),
    year = year(sample_week)
  ) |>
  summarise(flux = sum(mean, na.rm = T) * 7, .groups = "drop_last") |>
  summarise(AVG = mean(flux), SD = sd(flux))
message(
  "On average ",
  tot_flux |> filter(trophic_level == "zooplankton") |> pull(AVG) |> round(1),
  " kJ/m2 (SD = ",
  tot_flux |> filter(trophic_level == "zooplankton") |> pull(SD) |> round(1),
  ") are transfered from phytoplankton to zooplankton over the productive season"
)
percent_tot <- timeseries_fluxes |>
  add_season() |>
  filter(season != "Winter") |>
  group_by(
    trophic_level = ifelse(
      predator %in% node_data$node_name[node_data$trophic_level == 2],
      "zooplankton",
      "fish"
    ),
    year = year(sample_week)
  ) |>
  summarise(flux = sum(mean, na.rm = T) * 7, .groups = "drop") |>
  pivot_wider(values_from = flux, names_from = trophic_level) |>
  mutate(percent = 100 * fish / zooplankton)
message(
  "Between ",
  percent_tot |> pull(percent) |> min() |> round(1),
  " and ",
  percent_tot |> pull(percent) |> max() |> round(1),
  "% reached the three fish nodes"
)

# Fig S2 - 3: Impact on forage ratios on the fluxes ----
message("Generating Sup. Fig. S2 & S3")
# Impact of forage ratios
# or how the forage ratios deviate from the null model
isoweek_null <- timeseries_null |>
  filter(year(sample_week) > 2007) |>
  group_by(predator, prey, iso_week = isoweek(sample_week), station) |>
  summarise(avg_fluxes = mean(flux, na.rm = T), .groups = "drop")
ggsave(
  filename = file.path("output", "figure", "SupFigS2.pdf"),
  plot = plot_flux_difference("phytoplankton"),
  height = 10,
  width = 8,
  dpi = 500
)
ggsave(
  filename = file.path("output", "figure", "SupFigS3.pdf"),
  plot = plot_flux_difference("zooplankton") +
    scale_y_continuous(breaks = seq(0, 1, 0.001)),
  height = 8,
  width = 8,
  dpi = 500
)
