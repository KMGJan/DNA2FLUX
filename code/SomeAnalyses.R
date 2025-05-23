#!/usr/bin/env Rscript

# Process the fluxes if not done yet ----
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
  system(paste("nohup Rscript", file.path("code", "ProcessFluxes.R")))
}
rm(processedFluxes, fluxPaths)
# Load libraries ----
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
  library(data.table)
  library(spaa)
})

# Import the data ----
readAndArrange <- function(file = file, arrange = TRUE) {
  df <- read_csv(file = file, show_col_types = FALSE)
  if (arrange) df |> arrange(sample_week) else df
}
timeseries_null <- readAndArrange(
  file = file.path("data", "analyses", "timeseries_null.csv")
)
timeseries_fluxes <- readAndArrange(
  file = file.path("data", "analyses", "timeseries_fluxes.csv")
)
timeseries_nodes <- readAndArrange(
  file = file.path("data", "analyses", "timeseries_nodes.csv")
)
isoweek_fluxes <- readAndArrange(
  file = file.path("data", "analyses", "isoweek_fluxes.csv"),
  arrange = FALSE
)
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
)
# Set the color scheme for plotting ----
color_mapping <- setNames(node_data$color, node_data$node_name)
se <- function(x, ...) {
  n <- length(x[!is.na(x)]) # calculate the length of the vector
  if (n > 2) {
    # only compute standard error for vector >= 2
    out <- sd(x, ...) / sqrt(n)
  } else {
    out <- NA
  }
  return(out)
}
# Foodweb in spring, summer and fall ----

# Assign season based on ISO week
add_season <- function(df) {
  # Spring: week 5-25, Summer: 26-37, Fall: 38-50, Winter the remaining weeks
  df |>
    mutate(
      iso_week = isoweek(sample_week),
      season = case_when(
        iso_week %in% 5:25 ~ "Spring",
        iso_week %in% 26:37 ~ "Summer",
        iso_week %in% 38:50 ~ "Fall",
        TRUE ~ "Winter"
      )
    )
}

# Ensure season is a factor with correct order
factor_season <- function(df) {
  df |>
    mutate(
      season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))
    )
}

# Filter by iso week and year, optionally summarise biomass
filter_my_df <- function(df, summarise = TRUE) {
  df <- df |>
    filter(isoweek(sample_week) %in% 2:51, year(sample_week) > 2007)
  if (summarise) {
    df |>
      group_by(node_name, season, station_name) |>
      summarise(avg = mean(biomass, na.rm = TRUE), .groups = "drop")
  } else {
    df
  }
}
# Biomass per season
biomass_per_season <- weekly_biomasses |>
  add_season() |>
  filter_my_df(summarise = TRUE) |>
  filter(station_name == "BY31 LANDSORTSDJ") |>
  left_join(node_data, by = join_by(node_name)) |>
  mutate(node_name = factor(node_name, levels = node_data$node_name)) |>
  factor_season()

# Summarise the fluxes per season
flux_season <- function(df) {
  df |>
    filter_my_df(summarise = FALSE) |>
    add_season() |>
    group_by(predator, prey, station, season) |>
    summarise(flux = mean(mean, na.rm = TRUE), .groups = "drop") |>
    group_by(
      TL = ifelse(
        predator %in% c("Clupea", "Sprattus", "Gasterosteus"),
        "fish",
        "zooplankton"
      ),
      season
    )
}

# Season graph
season_graph <- timeseries_fluxes |>
  flux_season() |>
  mutate(flux = (flux / sum(flux)) * 100) |>
  ungroup() |>
  as_tbl_graph() |>
  activate(nodes) |>
  left_join(rename(node_data, name = node_name), by = join_by(name)) |>
  activate(edges) |>
  factor_season()

# Season fluxes summary (for annotation)
season_fluxes <- timeseries_fluxes |>
  flux_season() |>
  summarise(flux = sum(flux), .groups = "drop") |>
  mutate(y = ifelse(TL == "fish", 2.5, 1.5), x = 7) |>
  factor_season()

# Plot
season_foodweb <-
  ggraph(
    graph = season_graph,
    layout = "manual",
    x = horizontal_position,
    y = trophic_level
  ) +
  # Add the link between all species relative to their contribution to the total fluxes between TL
  geom_edge_link(mapping = aes(edge_width = flux, col = flux)) + #, alpha = flux)) +
  scale_edge_color_gradient2(
    low = "#f0f0f0",
    mid = "#737373",
    high = "#000000",
    limits = c(0, 25),
    midpoint = 10,
    name = "Contribution fluxes [%]"
  ) +
  scale_edge_width(
    range = c(.5, 5),
    limits = c(0, 25),
    name = "Contribution fluxes [%]"
  ) +
  #scale_edge_alpha_continuous(guide = "none", range = c(0.15, 1)) +
  # Add the total fluxes between TL as a label
  geom_label(
    data = season_fluxes,
    mapping = aes(
      x = x,
      y = y,
      label = paste(round(flux, 3), "kJ/d/m2", sep = "\n")
    ),
    size = 3
  ) +
  # Add the node relative to their average biomass at each season, coloured by node_name
  geom_point(
    data = biomass_per_season,
    mapping = aes(
      x = horizontal_position,
      y = trophic_level,
      size = avg,
      fill = node_name
    ),
    shape = 21
  ) +
  scale_fill_manual(values = color_mapping) +
  scale_size_continuous(guide = "none") +
  # Facet by season
  facet_wrap(. ~ season) +
  # Fix the axis limits and change the label on the y-axis
  coord_cartesian(xlim = c(0.8, 13.2), ylim = c(0.8, 3.2)) +
  scale_y_continuous(
    breaks = 1:3,
    labels = c("Phytoplankton", "Zooplankton", "Fish")
  ) +
  # Fix the theme
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL)

ggsave(
  plot = season_foodweb,
  filename = file.path("output", "figure", "foodweb_season.pdf"),
  height = 7,
  width = 11,
  dpi = 500
)

# Relative flux by month ----
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
  ungroup() |>
  pivot_wider(names_from = prey, values_from = rel_flux, values_fill = 0) |>
  pivot_longer(5:23, names_to = "prey", values_to = "rel_flux") |>
  mutate(prey = factor(prey, levels = node_data$node_name)) |>
  filter(isoweek(sample_week) %in% 2:51)

tot_fluxes <- timeseries_fluxes |>
  group_by(sample_week, predator, station) |>
  summarise(tot_flux = sum(mean, na.rm = T), .groups = "drop") |>
  mutate(
    type = ifelse(
      predator %in% c("Clupea", "Gasterosteus", "Sprattus"),
      "fish",
      "zooplankton"
    ),
    predator = factor(predator, levels = node_data$node_name)
  ) |>
  filter(isoweek(sample_week) %in% 2:51)
total_flux <- tot_fluxes |>
  group_by(iso_week = isoweek(sample_week), predator) |>
  summarise(
    consum = mean(tot_flux),
    consum_min = consum - sd(tot_flux),
    consum_min = ifelse(consum_min < 0, 0, consum_min),
    consum_max = consum + sd(tot_flux),
    .groups = "drop"
  ) |>
  mutate(predator = factor(predator, levels = node_data$node_name))
p1 <- rel_contribution_fluxes |>
  group_by(iso_week = isoweek(sample_week), predator, prey) |>
  summarise(
    y = mean(rel_flux),
    ymin = y - sd(rel_flux),
    ymin = ifelse(ymin < 0, 0, ymin),
    ymax = y + sd(rel_flux),
    .groups = "drop"
  ) |>
  mutate(
    predator = factor(predator, levels = node_data$node_name),
    prey = factor(prey, levels = node_data$node_name)
  ) |>
  ggplot(aes(
    x = iso_week,
    y = y,
    ymin = ymin,
    ymax = ymax,
    col = prey,
    fill = prey
  )) +
  #geom_ribbon(alpha = .2)+
  geom_area(stat = "identity") +
  facet_grid(predator ~ .) +
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb, #c("Jan", "Mar", "May", "Jul", "Sep", "Nov"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, .5),
    expand = c(0, 0)
  ) +

  theme_bw() +
  labs(x = NULL, y = "Diet composition")

p2 <- total_flux |>
  ggplot(aes(x = iso_week, y = consum, ymin = consum_min, ymax = consum_max)) +
  geom_ribbon(alpha = .2) +
  geom_line() +
  facet_grid(predator ~ ., scales = "free") +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb, #c("Jan", "Mar", "May", "Jul", "Sep", "Nov"),
    expand = c(0, 0)
  ) +
  theme_bw() +
  labs(x = NULL, y = "Daily consumption \n [kJ/d/m2]")

plot_prop <- p1 + p2 + plot_layout(guides = "collect")
ggsave(
  plot = plot_prop,
  filename = file.path("output", "figure", "proportionOverTime.pdf"),
  height = 8.5,
  width = 8.5
)
# Ordination ----

plot_fluxes <- function(type, rel_data, tot_data, color_mapping) {
  scaling_factor <- if (type == "fish") 25 else 1
  y_line_expr <- expr(tot_flux * !!scaling_factor)
  y_axis_transform <- if (type == "fish") ~ . / 25 else ~.

  ggplot() +
    geom_area(
      data = filter(rel_data, type == !!type),
      mapping = aes(x = sample_week, y = rel_flux, fill = prey),
      stat = "identity",
      alpha = 0.7
    ) +
    geom_line(
      data = filter(tot_data, type == !!type),
      mapping = aes(x = sample_week, y = !!y_line_expr),
      col = "black",
      linewidth = 1.2
    ) +
    facet_grid(predator ~ .) +
    scale_fill_manual(values = color_mapping) +
    scale_y_continuous(
      sec.axis = sec_axis(
        y_axis_transform,
        name = "Total consumption \n [kJ/day/m2]"
      ),
      expand = c(0, 0)
    ) +
    scale_x_date(date_breaks = "2 years", expand = c(0, 0)) +
    theme_bw() +
    labs(y = "Relative contribution", x = NULL)
}
fish_fluxes <- plot_fluxes(
  "fish",
  rel_contribution_fluxes,
  tot_fluxes,
  color_mapping
)
zp_fluxes <- plot_fluxes(
  "zooplankton",
  rel_contribution_fluxes,
  tot_fluxes,
  color_mapping
)
timeseries_flux_plot <- fish_fluxes /
  zp_fluxes +
  plot_layout(guides = "collect", heights = c(3, 8))
ggsave(
  plot = timeseries_flux_plot,
  filename = file.path("output", "figure", "timeseries_flux_plot.pdf"),
  width = 10,
  height = 10,
  dpi = 500
)

mat_zp <-
  rel_contribution_fluxes |>
  filter(type == "fish", rel_flux > 0) |>
  select(sample_week, station, predator, prey, rel_flux) |>
  pivot_wider(names_from = prey, values_from = rel_flux, values_fill = 0)
df_pp <-
  timeseries_fluxes |>
  group_by(
    sample_week,
    station,
    predator = factor(predator, levels = node_data$node_name),
    prey = factor(prey, levels = node_data$node_name),
    type = ifelse(
      predator %in% c("Clupea", "Gasterosteus", "Sprattus"),
      "fish",
      "zooplankton"
    )
  ) |>
  summarise(flux = sum(mean, na.rm = T), .groups = "drop") |>
  filter(type == "fish") |>
  # Join with zooplankton data that contains the prey-specific contributions of primary producers
  left_join(
    rel_contribution_fluxes |>
      filter(type == "zooplankton", rel_flux > 0) |>
      # Convert from long to wide: each column will be a primary producer
      pivot_wider(names_from = prey, values_from = rel_flux, values_fill = 0) |>
      # Rename zooplankton predator as 'prey' to match the fish dataset
      # (zooplankton are prey in the fish-level trophic step)
      select(sample_week, prey = predator, station, 5:16),
    by = join_by(sample_week, prey, station)
  ) |>
  # Scale the fluxes with the relative fluxes of primary consumers
  # This traces the path: fish ← zooplankton ← primary producers
  mutate(across(7:18, ~ . * flux)) |>
  # Remove the original rel_flux column — it's already incorporated in the scaling
  select(-flux) |>
  # Pivot from wide to long: each row represents a primary producer's contribution
  pivot_longer(
    6:17,
    names_to = "primary_producer",
    values_to = "scaled_rel_flux"
  ) |>
  # Group by sample week, fish predator, primary producer, and station
  # Then average contributions across zooplankton intermediates
  group_by(sample_week, predator, primary_producer, station) |>
  summarise(
    rel_primary_production = sum(scaled_rel_flux) / n_distinct(prey),
    .groups = "drop"
  ) |>
  # Normalize within each fish predator group (per week & station)
  group_by(sample_week, predator, station) |>
  mutate(
    rel_primary_production = rel_primary_production /
      sum(rel_primary_production)
  ) |>
  ungroup()

# PCoA + metadata helper
process_pcoa <- function(pcoa_result, metadata, id_cols = 1:3) {
  coords <- as.data.frame(pcoa_result$vectors[, 1:2])
  names(coords) <- c("Axis1", "Axis2")

  coords |>
    bind_cols(metadata[, id_cols]) |>
    mutate(
      iso_week = isoweek(sample_week),
      season = case_when(
        iso_week %in% 1:10 ~ "Winter",
        iso_week %in% 11:22 ~ "Spring",
        iso_week %in% 23:35 ~ "Summer",
        iso_week %in% 36:48 ~ "Autumn",
        iso_week %in% 49:53 ~ "Winter",
        TRUE ~ NA_character_
      ),
      month = month(sample_week),
      month_abb = factor(month.abb[month], levels = month.abb)
    )
}

# Prepare phytoplankton data (needs pivoting first)
plot_pcoa <- function(df, eig, envfit, title) {
  p1 <- ggplot(df, aes(x = Axis1, y = Axis2)) +

    geom_point(shape = 21, size = 1.5, mapping = aes(fill = predator)) +
    facet_wrap(~month_abb) +
    scale_fill_manual(values = color_mapping) +
    #coord_fixed() +
    theme_bw() +

    labs(
      title = title,
      x = paste0("Axis 1 (", round(100 * eig[1], 1), "%)"),
      y = paste0("Axis 2 (", round(100 * eig[2], 1), "%)")
    )
  p2 <- ggplot(envfit, mapping = aes(x = Axis1, y = Axis2)) +
    geom_segment(
      mapping = aes(x = Axis1, y = Axis2, yend = 0, xend = 0),
      arrow = arrow(length = unit(2, 'mm'), ends = "first")
    ) +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_x_continuous(limits = c(-1, 1)) +
    ggrepel::geom_text_repel(
      mapping = aes(x = Axis1, y = Axis2, label = prey),
      size = 2
    ) +
    coord_fixed() +
    theme_bw() +
    labs(
      x = paste0("Axis 1 (", round(100 * eig[1], 1), "%)"),
      y = paste0("Axis 2 (", round(100 * eig[2], 1), "%)")
    )
  p1 +
    p2 +
    plot_layout(guides = "collect", width = c(3, 1), axis_titles = "collect")
}
# For phytoplankton
mat_pp <- df_pp |>
  pivot_wider(
    names_from = primary_producer,
    values_from = rel_primary_production,
    values_fill = 0
  ) |>
  na.omit()
unique(isoweek(mat_pp$sample_week))
bray_dist_pp <- vegdist(mat_pp[, 4:15], method = "bray")
pcoa_result_pp <- pcoa(bray_dist_pp)
site_scores_pp <- pcoa_result_pp$vectors
envfit_pp <- as.data.frame(scores(
  envfit(site_scores_pp, mat_pp[4:15], permutations = 999),
  display = c("vectors")
)) |>
  rownames_to_column(var = "prey") |>
  rename("Axis1" = Axis.1, "Axis2" = Axis.2)
eig_pp <- pcoa_result_pp$values$Relative_eig
pcoa_df_pp <- process_pcoa(pcoa_result_pp, metadata = mat_zp)
pcoa_pp <- plot_pcoa(
  pcoa_df_pp,
  eig_pp,
  envfit_pp,
  "PCoA of Energy Fluxes based on relative contribution of phytoplankton to fish"
)

# For zooplankton
bray_dist_zp <- vegdist(mat_zp[, 4:10], method = "bray")
pcoa_result_zp <- pcoa(bray_dist_zp)
site_scores_zp <- pcoa_result_zp$vectors
envfit_zp <- as.data.frame(scores(
  envfit(site_scores_zp, mat_zp[4:10], permutations = 999),
  display = c("vectors")
)) |>
  rownames_to_column(var = "prey") |>
  rename("Axis1" = Axis.1, "Axis2" = Axis.2)
eig_zp <- pcoa_result_zp$values$Relative_eig
pcoa_df_zp <- process_pcoa(pcoa_result_zp, metadata = mat_zp)
pcoa_zp <- plot_pcoa(
  pcoa_df_zp,
  eig_zp,
  envfit_zp,
  "PCoA of Energy Fluxes based on relative contribution of phytoplankton to fish"
)
pcoa_plot <- pcoa_zp / pcoa_pp
ggsave(
  plot = pcoa_plot,
  filename = file.path("output", "figure", "pcoa.pdf"),
  dpi = 500,
  width = 10,
  height = 10
)

## PCOA gif ----
df <- pcoa_df_zp |>
  mutate(type = "zooplankton") |>
  bind_rows(
    pcoa_df_pp |>
      mutate(type = "phytoplankton")
  )
envfit <- envfit_zp |>
  mutate(type = "zooplankton") |>
  bind_rows(
    envfit_pp |>
      mutate(type = "phytoplankton")
  )

dir.create(file.path("output", "figure", "gif_frames"), recursive = TRUE)
for (week in 2:51) {
  month_abb <- df |> filter(iso_week == week) |> pull(month_abb) |> unique()
  if (length(month_abb) > 1) month_abb <- month_abb[1]
  plot <-
    ggplot(mapping = aes(x = Axis1, y = Axis2)) +
    geom_segment(
      data = envfit,
      mapping = aes(x = Axis1, y = Axis2, yend = 0, xend = 0),
      arrow = arrow(length = unit(2, 'mm'), ends = "first")
    ) +
    ggrepel::geom_text_repel(
      data = envfit,
      mapping = aes(x = Axis1, y = Axis2, label = prey),
      size = 2,
      seed = 100
    ) +
    stat_ellipse(
      data = df |> filter(iso_week == week),
      aes(col = predator, group = predator)
    ) +
    geom_point(
      data = df |> filter(iso_week == week),
      aes(fill = predator),
      shape = 21,
      size = 3
    ) +

    # geom_point(data = df |> filter(iso_week == week-1), alpha = .4, size = 3,
    #            aes(fill = predator), shape = 21)+
    # geom_point(data = df |> filter(iso_week == week-2), alpha = .1, size = 3,
    #            aes(fill = predator), shape = 21)+

    facet_grid(. ~ type) +
    coord_fixed(xlim = c(-.7, 1), ylim = c(-1, .7)) +
    theme_bw() +
    scale_fill_manual(values = color_mapping) +
    scale_color_manual(values = color_mapping) +
    annotate(
      "label",
      x = -.5,
      y = .7,
      label = paste("Week:", week, "\n Month:", month_abb),
      size = 3,
      label.size = 0.5,
      label.r = unit(0.15, "lines"),
      fill = "#ffffffcc",
      color = "black",
      fontface = "bold"
    )
  ggsave(
    sprintf("./output/figure/gif_frames/frame_%03d.png", week),
    plot = plot,
    width = 10,
    height = 5,
    dpi = 300
  )
}

png_files <- list.files(
  "./output/figure/gif_frames",
  pattern = "*.png",
  full.names = TRUE
)
gifski::gifski(
  png_files,
  gif_file = "./output/figure/PCOA.gif",
  width = 3000,
  height = 1500,
  delay = 0.2
)

unlink(file.path("output", "figure", "gif_frames"), recursive = TRUE)

# Permanovas ----
# This step takes some time about 5 min so it is better to save the output to come back to it later
#
if (file.exists(file.path("output", "table", "permanova_zp.rds"))) {
  permanova_zp <- read_rds(
    file = file.path("output", "table", "permanova_zp.rds")
  )
} else {
  set.seed(100)
  permanova_zp <- adonis2(
    formula = bray_dist_zp ~ predator + factor(month),
    data = mutate(mat_zp, month = month(sample_week)),
    permutations = 999,
    #strata = mat_zp$month,
    parallel = 8,
    by = "margin"
  )
  permanova_zp |>
    write_rds(file = file.path("output", "table", "permanova_zp.rds"))
}

if (file.exists(file.path("output", "table", "permanova_pp.rds"))) {
  permanova_zp <- read_rds(
    file = file.path("output", "table", "permanova_pp.rds")
  )
} else {
  set.seed(100)
  permanova_pp <- adonis2(
    formula = bray_dist_pp ~ predator + factor(month),
    data = mutate(mat_pp, month = month(sample_week)),
    permutations = 999,
    #strata = mat_pp$month,
    parallel = 8,
    by = "margin"
  )
  permanova_pp |>
    write_rds(file = file.path("output", "table", "permanova_pp.rds"))
}


p1 <- pcoa_df_zp |>
  pivot_longer(1:2, values_to = "PCOA_scores", names_to = "axis") |>
  ggplot(aes(
    x = sample_week,
    y = PCOA_scores,
    col = predator,
    group = interaction(predator, axis)
  )) +
  geom_line() +
  facet_grid(axis ~ .) +
  scale_color_manual(values = color_mapping) +
  theme_bw() +
  labs(title = "Based on zooplankton")
p2 <- pcoa_df_pp |>
  pivot_longer(1:2, values_to = "PCOA_scores", names_to = "axis") |>
  ggplot(aes(
    x = sample_week,
    y = PCOA_scores,
    col = predator,
    group = interaction(predator, axis)
  )) +
  geom_line() +
  facet_grid(axis ~ .) +
  scale_color_manual(values = color_mapping) +
  theme_bw() +
  labs(title = "Based on phytoplankton")

p <- p1 / p2 + plot_layout(guides = "collect", axes = "collect")
ggsave(
  plot = p,
  filename = file.path("output", "figure", "PCOA_timeseries.pdf"),
  width = 7,
  height = 5
)
p1 <- pcoa_df_zp |>
  pivot_longer(1:2, values_to = "PCOA_scores", names_to = "axis") |>
  group_by(iso_week, axis, predator) |>
  summarise(
    y = mean(PCOA_scores),
    ymin = y - sd(PCOA_scores),
    ymax = y + sd(PCOA_scores),
    .groups = "drop"
  ) |>
  ggplot(aes(
    x = iso_week,
    y = y,
    ymin = ymin,
    ymax = ymax,
    fill = predator,
    col = predator
  )) +
  geom_ribbon(alpha = .2) +
  geom_line() +
  facet_grid(axis ~ .) +
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  theme_bw() +
  labs(title = "Based on zooplankton", x = NULL, y = NULL) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  )
p2 <- pcoa_df_pp |>
  pivot_longer(1:2, values_to = "PCOA_scores", names_to = "axis") |>
  group_by(iso_week, axis, predator) |>
  summarise(
    y = mean(PCOA_scores),
    ymin = y - sd(PCOA_scores),
    ymax = y + sd(PCOA_scores),
    .groups = "drop"
  ) |>
  ggplot(aes(
    x = iso_week,
    y = y,
    ymin = ymin,
    ymax = ymax,
    fill = predator,
    col = predator
  )) +
  geom_ribbon(alpha = .2) +
  geom_line() +
  facet_grid(axis ~ .) +
  scale_color_manual(values = color_mapping) +
  scale_fill_manual(values = color_mapping) +
  theme_bw() +
  labs(title = "Based on phytoplankton", x = NULL, y = NULL) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  )
p <- p1 / p2 + plot_layout(guides = "collect", axes = "collect")
ggsave(
  plot = p,
  filename = file.path("output", "figure", "PCOA_timeseries_summary.pdf"),
  width = 5.5,
  height = 5
)
standardized_flux <- total_flux |>
  group_by(predator) |>
  mutate(z = (consum - mean(consum)) / sd(consum)) |>
  ungroup()

centroid_PCOA <-
  pcoa_df_pp |>
  mutate(type = "phytoplankton") |>
  bind_rows(pcoa_df_zp |> mutate(type = "zooplankton")) |>
  group_by(iso_week, predator, season, type) |>
  summarise(
    avg1 = mean(Axis1, na.rm = TRUE),
    avg2 = mean(Axis2, na.rm = TRUE),
    sd1 = sd(Axis1, na.rm = TRUE),
    sd2 = sd(Axis2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(standardized_flux, by = join_by(iso_week, predator)) |>
  ggplot() +
  geom_segment(
    data = envfit,
    mapping = aes(x = Axis1, y = Axis2, yend = 0, xend = 0),
    arrow = arrow(length = unit(2, 'mm'), ends = "first")
  ) +
  ggrepel::geom_text_repel(
    data = envfit,
    mapping = aes(x = Axis1, y = Axis2, label = prey),
    size = 2,
    seed = 100
  ) +
  geom_path(aes(x = avg1, y = avg2, col = predator, linewidth = z)) + #, linewidth = 1.5) +
  geom_point(
    aes(
      x = avg1,
      y = avg2,
      fill = iso_week,
      col = predator,
      group = predator,
      size = z
    ),
    #size = 3,
    shape = 21
  ) +
  facet_grid(. ~ type) +
  coord_fixed() +
  scale_color_manual(values = color_mapping) +
  scale_fill_gradientn(
    colors = c("#51291E", "#649A47", "#FEEA00", "#BCABAE", "#51291E"),
    values = scales::rescale(c(2, 5, 25, 37, 50)), # Rescales to [0,1]
    limits = c(1, 52), # Optional: match full range of iso_week
    breaks = c(12, 24, 36, 48),
    labels = c("Mar", "Jun", "Sep", "Dec")
  ) +
  theme_bw() +
  scale_size_continuous(
    name = "consumption \n z-score",
    limits = c(-1.5, 1.6),
    range = c(1, 5)
  ) +
  scale_linewidth_continuous(
    guide = "none",
    limits = c(-1.5, 1.6),
    range = c(1, 5)
  )


ggsave(
  plot = centroid_PCOA,
  filename = file.path("output", "figure", "centroid_PCOA.pdf"),
  dpi = 500,
  width = 10,
  height = 5
)

# Diet overlap -----
interspecific_overlap <- function(mat, cols, method = "schoener") {
  # Create numeric matrix and row names
  mat2 <- mat[, cols] |> as.matrix()
  rownames(mat2) <- with(mat, paste(station, predator, sample_week, sep = "_"))

  # Get Schoener's D overlap matrix and long format
  niche_df <- niche.overlap(t(mat2), method = method) |>
    dist2list() |>
    as.data.frame() |>
    mutate(value = ifelse(col == row, NA, as.numeric(value))) |>
    separate(
      col,
      into = c("station.x", "predator.x", "sample_week.x"),
      sep = "_"
    ) |>
    separate(
      row,
      into = c("station.y", "predator.y", "sample_week.y"),
      sep = "_"
    ) |>
    mutate(
      x = predator.x,
      y = predator.y,
      xy_min = pmin(x, y),
      xy_max = pmax(x, y)
    ) |>
    filter(x == xy_min, y == xy_max) |>
    select(-xy_min, -xy_max) |>
    na.omit()
}

plot_interspecific_overlap <- function(data) {
  data |>
    mutate(interaction = paste(x, y, sep = "_")) |>
    group_by(iso_week = isoweek(sample_week.x), interaction) |>
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
    geom_ribbon(alpha = .2) +
    geom_line() +
    theme_bw() +
    #scale_y_continuous(limits = c(0, 1)) +
    labs(x = NULL, y = "Schoener's D") +
    scale_x_continuous(
      breaks = seq(1, 52.1775, 4.348125),
      labels = month.abb,
      expand = c(0, 0)
    )
}
plot_zp <- mat_zp |>
  mutate(
    group_id = paste(year(sample_week), isoweek(sample_week), sep = "_")
  ) |>
  group_by(group_id) |>
  group_split() |>
  map_dfr(
    ~ interspecific_overlap(.x, 4:10)
  ) |>
  plot_interspecific_overlap()
plot_pp <- mat_pp |>
  mutate(
    group_id = paste(year(sample_week), isoweek(sample_week), sep = "_")
  ) |>
  group_by(group_id) |>
  group_split() |>
  map_dfr(
    ~ interspecific_overlap(.x, 4:15)
  ) |>
  plot_interspecific_overlap()
plot_schoeners <- plot_pp /
  plot_zp +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect")
ggsave(
  plot = plot_schoeners,
  filename = file.path("output", "figure", "overlap_across_weeks.pdf"),
  dpi = 500,
  width = 8,
  height = 5
)

## between week overlap ----
compute_group_overlap <- function(df, prey_cols) {
  if (nrow(df) < 2) return(NULL) # Skip groups with 1 row

  mat2 <- as.matrix(df[prey_cols])
  rownames(mat2) <- paste(df$station, df$predator, df$sample_week, sep = "_")

  niche.overlap(t(mat2), method = "schoener") |>
    dist2list() |>
    as.data.frame() |>
    mutate(value = ifelse(col == row, NA, as.numeric(value))) |>
    separate(
      col,
      into = c("station.x", "predator.x", "sample_week.x"),
      sep = "_"
    ) |>
    separate(
      row,
      into = c("station.y", "predator.y", "sample_week.y"),
      sep = "_"
    ) |>
    mutate(
      sample_week.x = as.Date(sample_week.x),
      sample_week.y = as.Date(sample_week.y)
    )
}

plot_month_overlap <- function(df, title) {
  df |>
    mutate(
      week.x = isoweek(sample_week.x),
      week.y = isoweek(sample_week.y),
      year = lubridate::year(sample_week.x),
      predator = predator.x,
      x = week.x,
      y = week.y,
      xy_min = pmax(x, y),
      xy_max = pmin(x, y)
    ) |>
    filter(x == xy_min, y == xy_max) |>
    group_by(week.y, week.x, predator, year) |>
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") |>
    group_by(week.x, week.y, predator) |>
    summarise(
      avg = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    na.omit() |>
    ggplot(aes(
      x = week.y,
      y = factor(week.x, levels = sort(unique(week.x), decreasing = TRUE)),
      fill = avg
    )) +
    geom_tile() +
    facet_grid(. ~ predator) +
    coord_fixed() +
    scale_fill_gradient2(
      mid = "#fdbb84",
      high = "black",
      low = "white",
      midpoint = .5,
      limits = c(0, 1),
      name = "Schoener's D"
    ) +
    theme_bw() +
    scale_x_continuous(
      breaks = seq(1, 52.1775, 4.348125),
      labels = month.abb,
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      breaks = round(seq(1, 52.1775, 4.348125)),
      labels = month.abb,
      expand = c(0, 0)
    ) +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = NULL, title = title)
}

overlap_zp <- mat_zp |>
  mutate(group_id = paste(predator, year(sample_week), sep = "_")) |>
  group_by(group_id) |>
  group_split() |>
  map_dfr(
    ~ compute_group_overlap(.x, 4:10)
  ) |>
  plot_month_overlap(title = "zooplankton")
overlap_pp <- mat_pp |>
  mutate(group_id = paste(predator, year(sample_week), sep = "_")) |>
  group_by(group_id) |>
  group_split() |>
  map_dfr(
    ~ compute_group_overlap(.x, 4:15)
  ) |>
  plot_month_overlap(title = "phytoplankton")
plot_overlap_over_month <- overlap_pp /
  overlap_zp +
  plot_layout(guides = "collect", axes = "collect")
ggsave(
  plot = plot_overlap_over_month,
  filename = file.path("output", "figure", "overlap_over_month.pdf"),
  width = 10,
  height = 7,
  dpi = 500
)
# Impact of forage ratios ----
# or how the forage ratios deviate from the null model
isoweek_null <- timeseries_null |>
  filter(year(sample_week) > 2007) |>
  group_by(predator, prey, iso_week = isoweek(sample_week), station) |>
  summarise(avg_fluxes = mean(flux, na.rm = T), .groups = "drop")

# Reusable plot function
plot_flux_difference <- function(trophic_level_filter) {
  isoweek_fluxes |>
    left_join(
      isoweek_null |> rename("null" = avg_fluxes),
      by = c("iso_week", "station", "predator", "prey")
    ) |>
    mutate(
      predator = factor(predator, levels = node_data$node_name),
      prey = factor(prey, levels = node_data$node_name),
      trophic_level = case_when(
        prey %in% node_data$node_name[node_data$trophic_level == 2] ~
          "zooplankton",
        prey %in% node_data$node_name[node_data$trophic_level == 1] ~
          "phytoplankton"
      )
    ) |>
    filter(trophic_level == trophic_level_filter, iso_week %in% 2:51) |>
    na.omit() |>
    ggplot(aes(x = iso_week)) +

    geom_ribbon(
      aes(y = mean, ymin = lower, ymax = upper),
      alpha = .2,
      linewidth = .2,
      fill = "black",
      col = "black"
    ) +

    geom_line(aes(y = mean, color = "Selectivity")) +
    geom_point(aes(y = null, color = "No selectivity"), size = .5) +
    geom_line(aes(y = mean, color = "Selectivity")) +

    facet_grid(
      prey ~ predator,
      scales = "free_y"
    ) +

    scale_color_manual(
      name = NULL,
      values = c("Selectivity" = "black", "No selectivity" = "#ff7f00")
    ) +

    scale_x_continuous(
      breaks = seq(1, 52.1775, 4.348125 * 2),
      labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"),
      expand = c(0, 0)
    ) +

    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    labs(
      y = "Fluxes \n [kJ/day/m2]",
      x = NULL
    )
}
phytoplankton_diff <- plot_flux_difference("phytoplankton")
zooplankton_diff <- plot_flux_difference("zooplankton")
seasonal_difference <- zooplankton_diff +
  phytoplankton_diff +
  plot_layout(width = c(3, 8), guides = "collect", axis_titles = "collect")
ggsave(
  filename = file.path("output", "figure", "seasonal_difference.pdf"),
  plot = seasonal_difference,
  height = 9,
  width = 15,
  dpi = 500
)
## Over the entire timeseries ----
tot_null <- timeseries_null |>
  filter(year(sample_week) > 2007) |>
  group_by(station, predator, prey) |>
  summarise(flux = mean(flux, na.rm = T), .groups = "drop")
flux_diff <- station_fluxes |>
  left_join(tot_null, by = c("station", "predator", "prey")) |>
  mutate(
    sig = ifelse(
      (lower - flux > 0) | (upper - flux < 0),
      "Significant",
      "Not significant"
    ),
    prey = factor(prey, levels = node_data$node_name) |> fct_rev(),
    predator = factor(predator, levels = node_data$node_name),
    trophic_level = case_when(
      prey %in% node_data$node_name[node_data$trophic_level == 2] ~
        "zooplankton",
      prey %in% node_data$node_name[node_data$trophic_level == 1] ~
        "phytoplankton"
    )
  ) |>
  group_by(station, predator, trophic_level, sig) |>
  mutate(
    rel_mean = (mean - flux), # / flux * 100,
    rel_upper = (upper - flux), # / flux * 100,
    rel_lower = (lower - flux), # / flux * 100,
    rel_null = 0,
    prey = factor(prey, levels = node_data$node_name),
    predator = factor(predator, levels = node_data$node_name),
    predator = fct_rev(predator)
  ) |>
  ungroup() |>
  filter(!is.na(flux))

# Plotting function
plot_flux_diff <- function(data, level) {
  data |>
    filter(trophic_level == level) |>
    ggplot(aes(
      x = prey,
      y = rel_mean,
      ymax = rel_upper,
      ymin = rel_lower,
      fill = predator,
      alpha = sig
    )) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +

    geom_bar(
      stat = "identity",
      #shape = 21,
      #size = 2,
      position = position_dodge2(width = 0.4, preserve = "single")
    ) +
    geom_errorbar(
      position = position_dodge2(width = 0.4, preserve = "single")
    ) +
    #facet_grid(. ~ predator, scales = "fixed") +
    scale_fill_manual(values = color_mapping) +
    scale_alpha_manual(values = c(0.4, 1)) +
    labs(y = "Difference from null model [kJ/day/m2]", x = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

plot_zoop <-
  plot_flux_diff(flux_diff, "zooplankton") +
  scale_y_continuous(limits = c(-0.0035, 0.0035))
plot_phyto <- plot_flux_diff(flux_diff, "phytoplankton") +
  scale_y_continuous(limits = c(-0.02, 0.02))
plot_diff <- plot_zoop +
  plot_phyto +
  plot_layout(width = c(3, 7), guides = "collect", axes = "collect")
ggsave(
  plot = plot_diff,
  filename = file.path("output", "figure", "annual_differences.pdf"),
  width = 10,
  height = 5,
  dpi = 500
)

# Predation pressure ----
# Join and calculate predation pressure
timeseries_pressure <-
  timeseries_nodes |>
  filter(!name %in% c("Clupea", "Gasterosteus", "Sprattus")) |>
  rename("prey" = name) |>
  select(prey, sample_week, biomass) |>
  right_join(
    timeseries_fluxes,
    by = c("prey", "sample_week"),
    relationship = "many-to-many"
  ) |>

  group_by(
    prey,
    predator,
    sample_week,
    station,
    iso_week = isoweek(sample_week),
    year,
    biomass
  ) |>
  summarise(pressure_sum = sum(mean, na.rm = TRUE), .groups = "drop_last") |>
  summarise(
    pressure = sum(pressure_sum, na.rm = TRUE) / sum(biomass, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(iso_week %in% 2:51) |>
  mutate(
    trophic_level = case_when(
      prey %in% node_data$node_name[node_data$trophic_level == 2] ~
        "zooplankton",
      prey %in% node_data$node_name[node_data$trophic_level == 1] ~
        "phytoplankton"
    )
  ) |>
  pivot_wider(names_from = prey, values_from = pressure, values_fill = 0) |>
  pivot_longer(7:25, names_to = "prey", values_to = "pressure") |>
  group_by(predator, prey) |>
  filter(sum(pressure) > 0) |>
  ungroup()

timeseries_pressure |>
  group_by(iso_week, predator, prey) |>
  summarise(pressure = mean(pressure), .groups = "drop") |>
  mutate(
    prey = factor(prey, levels = node_data$node_name),
    predator = factor(predator, levels = node_data$node_name)
  ) |>
  ggplot(aes(x = iso_week, y = pressure, fill = prey)) +
  geom_area() +
  facet_grid(predator ~ ., scales = "free") +
  scale_fill_manual(values = c("Total" = "black", color_mapping)) +
  labs(y = "Predation pressure [kJ/day/g]", x = NULL) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  theme_bw()

timeseries_pressure |>
  group_by(iso_week, predator, prey) |>
  summarise(pressure = mean(pressure), .groups = "drop") |>
  mutate(
    prey = factor(prey, levels = node_data$node_name),
    predator = factor(predator, levels = node_data$node_name)
  ) |>
  ggplot(aes(x = iso_week, y = pressure, fill = predator)) +
  geom_area() +
  facet_grid(prey ~ ., scales = "free") +
  scale_fill_manual(values = c("Total" = "black", color_mapping)) +
  labs(y = "Predation pressure [kJ/day/g]", x = NULL) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  theme_bw()
timeseries_nodes |>
  filter(!name %in% c("Clupea", "Gasterosteus", "Sprattus")) |>
  rename("prey" = name) |>
  select(prey, sample_week, biomass) |>
  right_join(
    timeseries_fluxes,
    by = c("prey", "sample_week"),
    relationship = "many-to-many"
  ) |>

  group_by(
    predator,
    sample_week,
    station,
    iso_week = isoweek(sample_week),
    year,
    biomass
  ) |>
  summarise(pressure_sum = sum(mean, na.rm = TRUE), .groups = "drop_last") |>
  summarise(
    pressure = sum(pressure_sum, na.rm = TRUE) / sum(biomass, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(iso_week %in% 2:51) |>
  mutate(
    trophic_level = case_when(
      predator %in% node_data$node_name[node_data$trophic_level == 3] ~ "fish",
      predator %in% node_data$node_name[node_data$trophic_level == 2] ~
        "zooplankton"
    )
  ) |>
  ggplot(aes(x = sample_week, y = pressure, col = predator)) +
  geom_line() +
  facet_wrap(trophic_level ~ predator, scales = "free")

pred_pressure_plot <- timeseries_nodes |>
  filter(!name %in% c("Clupea", "Gasterosteus", "Sprattus")) |>
  rename("prey" = name) |>
  select(prey, sample_week, biomass, trophic_level) |>
  right_join(
    timeseries_fluxes,
    by = c("prey", "sample_week"),
    relationship = "many-to-many"
  ) |>

  group_by(
    trophic_level = case_when(
      trophic_level == 3 ~ "fish",
      trophic_level == 1 ~ "phytoplankton",
      TRUE ~ "zooplankton"
    ),
    sample_week,
    station,
    iso_week = isoweek(sample_week),
    year,
    biomass
  ) |>
  summarise(pressure_sum = sum(mean, na.rm = TRUE), .groups = "drop_last") |>
  summarise(
    pressure = sum(pressure_sum, na.rm = TRUE) / sum(biomass, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    pressure = ifelse(trophic_level == "zooplankton", pressure * 2.5, pressure)
  ) |>
  filter(iso_week %in% 2:51) |>
  group_by(iso_week, trophic_level, station) |>
  summarise(
    median = median(pressure),
    min = quantile(pressure, .25),
    max = quantile(pressure, .75),
    .groups = "drop"
  ) |>
  ggplot(aes(
    x = iso_week,
    y = median,
    ymin = min,
    ymax = max,
    col = trophic_level,
    fill = trophic_level
  )) +
  geom_ribbon(alpha = .5) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("#1e713f", "#635255")) +
  scale_fill_manual(values = c("#1e713f", "#635255")) +
  theme_bw() +
  labs(x = NULL, y = "predation pressure on phytoplankton") +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . / 2.5, name = "predation pressure on zooplankton")
  )
temperature <- readAndArrange(
  file = file.path("data", "processed", "interpolation", "temperature.csv")
)
temp <- temperature |>
  mutate(temperature = temperature * 1.9) |>
  group_by(iso_week = isoweek(sample_week), station_name) |>
  summarise(
    median = median(temperature),
    min = quantile(temperature, .25),
    max = quantile(temperature, .75),
    .groups = "drop"
  ) |>
  filter(station_name == "BY31 LANDSORTSDJ", iso_week %in% 2:51)
biomass_plot <-
  weekly_biomasses |>
  left_join(node_data, by = join_by(node_name)) |>
  filter(isoweek(sample_week) %in% 2:51, type != "fish") |>
  group_by(sample_week, type, station_name) |>
  summarise(biomass = sum(biomass), .groups = "drop") |>
  group_by(iso_week = isoweek(sample_week), type, station_name) |>
  summarise(
    median = median(biomass),
    min = quantile(biomass, .25),
    max = quantile(biomass, .75),
    .groups = "drop"
  ) |>
  filter(station_name == "BY31 LANDSORTSDJ") |>
  rename("trophic_level" = type) |>
  ggplot(aes(x = iso_week, y = median, ymin = min, ymax = max)) +
  geom_ribbon(data = temp, fill = "#6A041D", col = "#6A041D", alpha = .5) +
  geom_line(data = temp, col = "#6A041D", linewidth = 1) +
  geom_ribbon(
    alpha = .5,
    mapping = aes(fill = trophic_level, col = trophic_level)
  ) +
  geom_line(linewidth = 1, mapping = aes(col = trophic_level)) +

  scale_color_manual(values = c("#1e713f", "#635255")) +
  scale_fill_manual(values = c("#1e713f", "#635255")) +
  theme_bw() +
  labs(x = NULL, y = "Biomass [g/m2]") +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  scale_y_continuous(sec.axis = sec_axis(~ . / 1.9, name = "Temperature [°C]"))
combined_plot <- pred_pressure_plot /
  biomass_plot +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect")
ggsave(
  plot = combined_plot,
  filename = file.path("output", "figure", "predation_pressure.pdf"),
  height = 5,
  width = 7,
  dpi = 500
)
# Annual fluxes ----
annual_plot <- timeseries_fluxes |>
  left_join(node_data, by = join_by(predator == "node_name")) |>
  group_by(sample_week, iso_week, trophic_level, station) |>
  summarise(flux = sum(mean, na.rm = T), .groups = "drop") |>
  mutate(
    interaction = ifelse(
      trophic_level == 2,
      "zooplankton consumption",
      "fish consumption"
    ),
    interaction = factor(
      interaction,
      levels = c("fish consumption", "zooplankton consumption")
    )
  ) |>
  filter(year(sample_week) > 2007) |>
  mutate(flux = ifelse(trophic_level == 3, flux * 4.3, flux)) |>
  group_by(interaction, year = year(sample_week), station) |>
  summarise(flux = mean(flux, na.rm = T) * 365, .groups = "drop") |>
  group_by(
    before_wave = case_when(year < 2011 ~ "before", TRUE ~ "after"),
    interaction
  ) |>
  mutate(avg = mean(flux)) |>
  ungroup() |>
  ggplot(aes(x = year, y = flux, fill = interaction, col = interaction)) +
  geom_smooth(method = "loess", formula = "y ~ x") +
  geom_point(shape = 21, size = 3, col = "black") +
  labs(y = "zooplankton annual fluxes \n [kJ/yr/m2]", x = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  scale_x_continuous(breaks = seq(2008, 2030, 2)) +
  scale_color_manual(values = c("#2F195F", "#635255")) +
  scale_fill_manual(values = c("#2F195F", "#635255")) +
  scale_y_continuous(
    breaks = seq(0, 125, 25),
    sec.axis = sec_axis(
      ~ . / 4.3,
      name = "fish annual consumption [kJ/yr/m2]",
      breaks = seq(0, 50, 5)
    )
  )
ggsave(
  plot = annual_plot,
  filename = file.path("output", "figure", "annual_fluxes.pdf"),
  width = 6.5,
  height = 4,
  dpi = 500
)
annual_efficiency <-
  timeseries_fluxes |>
  left_join(node_data, by = join_by(predator == "node_name")) |>
  group_by(sample_week, iso_week, trophic_level, station) |>
  summarise(flux = sum(mean, na.rm = T), .groups = "drop") |>
  mutate(
    interaction = ifelse(
      trophic_level == 2,
      "zooplankton consumption",
      "fish consumption"
    ),
    interaction = factor(
      interaction,
      levels = c("fish consumption", "zooplankton consumption")
    )
  ) |>
  filter(year(sample_week) > 2007) |>
  group_by(interaction, year = year(sample_week), station) |>
  summarise(flux = mean(flux, na.rm = T) * 365, .groups = "drop") |>
  pivot_wider(names_from = interaction, values_from = flux) |>
  mutate(efficiency = 100 * `fish consumption` / `zooplankton consumption`) |>
  ggplot(aes(x = year, y = efficiency)) +
  geom_point(shape = 21, size = 3, col = "black") +
  geom_smooth(method = "loess", formula = "y ~ x", col = "black") +
  theme_bw() +
  labs(x = NULL, y = "trophic efficiency [%]")
ggsave(
  plot = annual_efficiency,
  filename = file.path("output", "figure", "annual_efficiency.pdf"),
  width = 6.5,
  height = 3,
  dpi = 500
)
timeseries_fluxes |>
  left_join(node_data, by = join_by(predator == "node_name")) |>
  group_by(sample_week, iso_week, trophic_level, station) |>
  summarise(flux = sum(mean, na.rm = T), .groups = "drop") |>
  mutate(
    interaction = ifelse(
      trophic_level == 2,
      "zooplankton consumption",
      "fish consumption"
    ),
    interaction = factor(
      interaction,
      levels = c("fish consumption", "zooplankton consumption")
    )
  ) |>
  filter(year(sample_week) > 2007) |>
  group_by(interaction, year = year(sample_week), station) |>
  summarise(flux = mean(flux, na.rm = T) * 365, .groups = "drop") |>
  pivot_wider(names_from = interaction, values_from = flux) |>
  mutate(efficiency = 100 * `fish consumption` / `zooplankton consumption`) |>
  pull(efficiency) |>
  mean()
# Ratio between Trophic Levels ----
timeseries_fluxes |>
  left_join(node_data, by = c(predator = "node_name")) |>
  group_by(sample_week, type, station) |>
  summarise(flux = sum(mean, na.rm = T), .groups = "drop") |>
  pivot_wider(names_from = type, values_from = flux) |>
  mutate(ratio = fish / zooplankton) |>
  ggplot(aes(x = sample_week)) +
  geom_line(aes(y = fish), col = "#2F195F") +
  geom_line(aes(y = zooplankton), col = "#635255") +
  geom_line(aes(y = ratio)) +
  theme_bw()

iso_week_fluxes_with_ratio <-
  timeseries_fluxes |>
  left_join(node_data, by = c(predator = "node_name")) |>
  group_by(sample_week, type, station) |>
  summarise(flux = sum(mean, na.rm = T), .groups = "drop") |>
  pivot_wider(names_from = type, values_from = flux) |>
  mutate(ratio = fish / zooplankton) |>
  pivot_longer(3:5) |>
  mutate(value = ifelse(name == "fish", value * 20, value)) |>
  group_by(iso_week = isoweek(sample_week), station, name) |>
  summarise(
    avg = mean(value),
    min = avg - sd(value),
    max = avg + sd(value),
    .groups = "drop"
  ) |>
  filter(iso_week %in% 2:51) |>
  mutate(min = ifelse(min < 0, 0, min))

flux_plot <-
  ggplot(mapping = aes(x = iso_week, y = avg, ymin = min, ymax = max)) +
  geom_ribbon(
    data = iso_week_fluxes_with_ratio |> filter(name == "fish"),
    col = "#FFC759",
    fill = "#FFC759",
    alpha = .5
  ) +
  geom_ribbon(
    data = iso_week_fluxes_with_ratio |> filter(name == "zooplankton"),
    col = "#607196",
    fill = "#607196",
    alpha = .5
  ) +
  geom_line(
    data = iso_week_fluxes_with_ratio |> filter(name == "fish"),
    col = "#FFC759"
  ) +
  geom_line(
    data = iso_week_fluxes_with_ratio |> filter(name == "zooplankton"),
    col = "#607196"
  ) +
  theme_bw() +
  scale_y_continuous(
    name = "zooplankton consumption \n [kJ/d/m2]",
    sec.axis = sec_axis(~ . / 20, name = "fish consumption \n [kJ/d/m2]"),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  labs(x = NULL) +
  theme(
    axis.title.y.left = element_text(color = "#607196"),
    axis.title.y.right = element_text(color = "#FFC759")
  )

ratio_plot <-
  ggplot(mapping = aes(x = iso_week, y = avg, ymin = min, ymax = max)) +
  geom_ribbon(
    data = iso_week_fluxes_with_ratio |> filter(name == "ratio"),
    col = "#4E3822",
    fill = "#4E3822",
    alpha = .5
  ) +
  geom_line(
    data = iso_week_fluxes_with_ratio |> filter(name == "ratio"),
    col = "#4E3822"
  ) +
  theme_bw() +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "contribution of fish predation \n on zooplankton consumption"
  ) +
  scale_y_continuous(
    limits = c(0, .75),
    breaks = seq(0, 75, .25),
    expand = c(0, 0)
  ) +
  theme(axis.title.y.left = element_text(color = "#4E3822"))

fluxes_and_ratio <- flux_plot / ratio_plot + plot_layout(axes = "collect")

ggsave(
  plot = fluxes_and_ratio,
  filename = file.path("output", "figure", "fluxes_and_ratio.pdf"),
  width = 5,
  height = 4
)
## Species specific ratios ----
zp_specific <-
  timeseries_fluxes |>
  left_join(node_data, by = c(predator = "node_name")) |>
  group_by(
    sample_week,
    type,
    station,
    prey = ifelse(type == "zooplankton", "phytoplankton", prey)
  ) |>
  summarise(flux = sum(mean, na.rm = T), .groups = "drop")
zp_fluxes <- zp_specific |>
  filter(type == "zooplankton") |>
  rename("zooplankton_consumption" = flux) |>
  select(-prey) |>
  left_join(
    zp_specific |>
      filter(type != "zooplankton") |>
      select(-type),
    by = join_by(sample_week, station)
  ) |>
  select(sample_week, station, prey, flux, zooplankton_consumption)
zp_fluxes_and_ratio <-
  zp_fluxes |>
  mutate(ratio = flux / zooplankton_consumption) |>
  select(-c(flux, zooplankton_consumption)) |>
  pivot_wider(names_from = prey, values_from = ratio, values_fill = 0) |>
  pivot_longer(3:9, values_to = "ratio") |>
  left_join(
    zp_fluxes |>
      select(-zooplankton_consumption) |>
      pivot_wider(names_from = prey, values_from = flux, values_fill = 0) |>
      pivot_longer(3:9, values_to = "flux"),
    by = join_by(sample_week, station, name)
  ) |>
  mutate(flux = flux * 20) |>
  group_by(iso_week = isoweek(sample_week), station, name) |>
  summarise(
    ratio_avg = mean(ratio, na.rm = T),
    ratio_min = ratio_avg - sd(ratio, na.rm = T),
    ratio_max = ratio_avg + sd(ratio, na.rm = T),
    consum_avg = mean(flux, na.rm = T),
    consum_min = consum_avg - sd(flux, na.rm = T),
    consum_max = consum_avg + sd(flux, na.rm = T),
    .groups = "drop"
  ) |>
  filter(iso_week %in% 2:51) |>
  mutate(
    ratio_min = ifelse(ratio_min < 0, 0, ratio_min),
    consum_min = ifelse(consum_min < 0, 0, consum_min),
    name = factor(name, levels = node_data$node_name)
  ) |>
  ggplot(aes(x = iso_week)) +
  geom_ribbon(
    alpha = .5,
    col = "#4E3822",
    fill = "#4E3822",
    mapping = aes(y = ratio_avg, ymin = ratio_min, ymax = ratio_max)
  ) +
  geom_ribbon(
    alpha = .5,
    col = "#607196",
    fill = "#607196",
    mapping = aes(y = consum_avg, ymin = consum_min, ymax = consum_max)
  ) +
  geom_line(mapping = aes(y = ratio_avg), col = "#4E3822") +
  geom_line(col = "#607196", mapping = aes(y = consum_avg)) +
  facet_grid(name ~ .) +
  theme_bw() +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "contribution of fish predation \n on zooplankton consumption"
  ) +
  theme(
    axis.title.y.left = element_text(color = "#4E3822"),
    axis.title.y.right = element_text(color = "#607196")
  ) +
  scale_y_continuous(
    limits = c(0, .5),
    breaks = seq(0, 75, .2),
    expand = c(0, 0),
    sec.axis = sec_axis(
      ~ . / 20,
      name = "zooplankton consumption \n [kJ/d/m2]",
      breaks = seq(0, .02, .01)
    )
  )
ggsave(
  filename = file.path("output", "figure", "zp_flux_and_ratio.pdf"),
  plot = zp_fluxes_and_ratio,
  height = 7,
  width = 4.3
)
