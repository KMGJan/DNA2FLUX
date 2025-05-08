#!/usr/bin/env Rscript

# Process the fluxes if not done yet ----
processedFluxes <- c("annual_fluxes.csv", "as_tbl_graph_timeseries_fluxes.rds", "isoweek_fluxes.csv", "station_fluxes.csv", "timeseries_fluxes.csv","timeseries_nodes.csv","timeseries_null.csv")
fluxPaths <- file.path("data", "analyses", processedFluxes)
if (!all(file.exists(fluxPaths))) {
  system(paste("nohup Rscript", file.path("code", "ProcessFluxes.R")))
}
rm(processedFluxes, fluxPaths)
# Load libraries ----
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(patchwork))

# Import the data ----
timeseries_null <- read_csv(file = file.path("data", "analyses", "timeseries_null.csv"), show_col_types = FALSE)
timeseries_fluxes <- read_csv(file = file.path("data", "analyses", "timeseries_fluxes.csv"), show_col_types = FALSE)
#timeseries_nodes <- read_csv(file = file.path("data", "analyses", "timeseries_nodes_with_growth.csv"), show_col_types = FALSE)
isoweek_fluxes <- read_csv(file = file.path("data", "analyses", "isoweek_fluxes.csv"), show_col_types = FALSE)
annual_fluxes <- read_csv(file = file.path("data", "analyses", "annual_fluxes.csv"), show_col_types = FALSE)
station_fluxes <- read_csv(file = file.path("data", "analyses", "station_fluxes.csv"), show_col_types = FALSE)
node_data <- read_csv(file = file.path("data", "raw", "node_data.csv"), show_col_types = FALSE)
weekly_biomasses <- read_csv(file = file.path("data", "processed","interpolation" , "weekly_biomasses.csv"), show_col_types = FALSE)
# Set the color scheme for plotting ----
color_mapping = setNames(node_data$color, node_data$node_name)

# Foodweb in spring, summer and fall ----
# Spring bloom: week 5-22. 23 and 37 summer bloom. 38 to 50 fall bloom
library(tidygraph)
biomass_per_season  <- timeseries_nodes |> 
  mutate(iso_week = isoweek(sample_week),
         season = ifelse(iso_week %in% 5:25, "Spring", ifelse(iso_week %in% 26:37, "Summer ", ifelse(iso_week %in% 38:50, "Fall", "Winter")))) |> 
  group_by(name, season) |> 
  summarise(biomass = mean(biomass, na.rm = T), .groups = "drop")


season_graph <-
  timeseries_fluxes |> 
  filter(iso_week %in% 2:51) |> 
  mutate(month = month(sample_week),
         season = ifelse(iso_week %in% 12:25, "Spring", ifelse(iso_week %in% 26:38, "Summer", ifelse(iso_week %in% 39:51, "Fall", "Winter")))) |> 
  group_by(predator, prey,station,  season) |> 
  summarise(flux = mean(mean, na.rm = T), .groups = "drop") |> 
  group_by(TL = ifelse(predator %in% c("Clupea", "Sprattus", "Gasterosteus"), "fish", "zooplankton"),
    #predator,
    season) |> 
  mutate(flux = flux / sum(flux)) |> 
  ungroup() |> 
  as_tbl_graph() |> 
  activate(nodes) |> 
  left_join(node_data |> rename("name" = node_name), by = join_by(name)) |> 
  activate(edges) |> 
  mutate(season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter")))

season_fluxes <-
  timeseries_fluxes |> 
  filter(iso_week %in% 2:51) |> 
  mutate(month = month(sample_week),
         season = ifelse(iso_week %in% 12:25, "Spring", ifelse(iso_week %in% 26:38, "Summer", ifelse(iso_week %in% 39:51, "Fall", "Winter")))) |> 
  group_by(predator, prey,station,  season) |> 
  summarise(flux = mean(mean, na.rm = T), .groups = "drop") |> 
  group_by(TL = ifelse(predator %in% c("Clupea", "Sprattus", "Gasterosteus"), "fish", "zooplankton"),
           #predator,
           season) |> 
  summarise(flux = sum(flux), .groups = "drop") |> 
  mutate(y = ifelse(TL == "fish", 2.5, 1.5), x = ((1+13)/2)) |> 
  mutate(season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter")))

season_graph |> 
  ggraph(layout = "manual", x = horizontal_position, y = trophic_level) +
  geom_edge_link(aes(edge_width = flux, col = flux, alpha = flux),
                 arrow = arrow(length = unit(2, 'mm'), ends = "first")) +

  geom_node_label(aes(label = name), angle = 45, size = 3, nudge_y = -0.05, nudge_x = .2, hjust = 1) +
  geom_label(data = season_fluxes, mapping = aes(x = x, y=y, label = paste(round(flux, 3), "kJ/d/m2", sep = "\n")), size = 3) +
  coord_cartesian(xlim = c(0,14),
                  ylim = c(0,3.5)) +
  scale_y_continuous(breaks = 1:3) +
  facet_wrap(~season) +
  theme_bw() +
  labs(x = NULL,
       y = "Trophic level") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  scale_edge_color_gradientn(colours = c("black", "#A23C2A", "#136F63"),
                             values = scales::rescale(c(0, 0.1, 0.25), from = c(0, 0.25)),
                             limits = c(0, 0.25), name = "Contribution fluxes [%]")+
  scale_edge_width(range = c(.5, 5), limits = c(0, 0.25), name = "Contribution fluxes [%]") +
  scale_edge_alpha_continuous(guide = "none", range = c(0.15, 1))
ggsave("./output/figure/foodweb_season.pdf", height = 10, width = 12, dpi = 500)

# Ordination ----
rel_contribution_fluxes <-
  timeseries_fluxes |>
  group_by(sample_week, predator) |>
  reframe(sample_week = sample_week,
          station = station,
          prey = factor(prey, levels = node_data$node_name),
          predator =factor(predator, levels = node_data$node_name),
          type = ifelse(predator %in% c("Clupea", "Gasterosteus", "Sprattus"), "fish", "zooplankton"),
          rel_flux = mean / sum(mean, na.rm = T)) |> 
  ungroup() |> 
  pivot_wider(names_from = prey, values_from = rel_flux, values_fill = 0) |>
  pivot_longer(5:23, names_to = "prey", values_to = "rel_flux") |> 
  mutate(prey = factor(prey, levels = node_data$node_name))# |> filter(prey == "Peridiniales")

tot_fluxes <- timeseries_fluxes |>
  group_by(sample_week, predator, station) |>
  summarise(tot_flux = sum(mean, na.rm = T), .groups = "drop") |> 
  mutate(type = ifelse(predator %in% c("Clupea", "Gasterosteus", "Sprattus"), "fish", "zooplankton"),
         predator =factor(predator, levels = node_data$node_name)) |> 
  filter(isoweek(sample_week) %in% 2:50)

plot_fluxes <- function(type, rel_data, tot_data, color_mapping) {
  scaling_factor <- if (type == "fish") 25 else 1/4
  y_line_expr <- expr(tot_flux * !!scaling_factor)
  y_axis_transform <- if (type == "fish") ~ . / 25 else ~ . * 4
  
  ggplot() +
    geom_area(data = filter(rel_data, type == !!type),
              mapping = aes(x = sample_week, y = rel_flux, fill = prey),
              stat = "identity", alpha = 0.7) +
    geom_line(data = filter(tot_data, type == !!type),
              mapping = aes(x = sample_week, y = !!y_line_expr ),
              col = "black", linewidth = 1.2) +
    facet_grid(predator ~ .) +
    scale_fill_manual(values = color_mapping) +
    scale_y_continuous(sec.axis = sec_axis(y_axis_transform,
                                           name = "Total incoming fluxes \n [kJ/day/m2]"),
                       expand = c(0,0)) +
    scale_x_date(date_breaks = "2 years", expand = c(0,0)) +
    theme_bw()+
    labs(y = "Relative contribution",
         x = NULL)
}
fish_fluxes <- plot_fluxes("fish", rel_contribution_fluxes, tot_fluxes, color_mapping)
zp_fluxes <- plot_fluxes("zooplankton", rel_contribution_fluxes, tot_fluxes, color_mapping)
fish_fluxes / zp_fluxes + plot_layout(guides = "collect", heights = c(3, 8))


mat_zp <- 
  rel_contribution_fluxes |>
  filter(type == "fish",
         rel_flux > 0) |>
  select(sample_week, station, predator, prey, rel_flux) |> 
  pivot_wider(names_from = prey, values_from = rel_flux, values_fill = 0)
df_pp <-
  timeseries_fluxes |>
  group_by(sample_week,
           station,
           predator = factor(predator, levels = node_data$node_name),
           prey = factor(prey, levels = node_data$node_name),
           type = ifelse(predator %in% c("Clupea", "Gasterosteus", "Sprattus"), "fish", "zooplankton")) |> 
  summarise(flux = sum(mean, na.rm = T), .groups = "drop") |>
  filter(type == "fish") |> 
  # Join with zooplankton data that contains the prey-specific contributions of primary producers
  left_join(
    rel_contribution_fluxes |> 
      filter(type == "zooplankton",
             rel_flux > 0) |> 
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
  pivot_longer(6:17, names_to = "primary_producer", values_to = "scaled_rel_flux") |> 
  # Group by sample week, fish predator, primary producer, and station
  # Then average contributions across zooplankton intermediates
  group_by(sample_week, predator, primary_producer, station) |> 
  summarise(rel_primary_production = sum(scaled_rel_flux) / n_distinct(prey), .groups = "drop") |> 
  # Normalize within each fish predator group (per week & station)
  group_by(sample_week, predator, station) |> 
  mutate(rel_primary_production = rel_primary_production/sum(rel_primary_production)) |> 
  ungroup()

df_pp |> 
  ggplot(aes(x= sample_week, y = rel_primary_production, fill = primary_producer))+
  geom_area(stat = "identity")+
  facet_grid(predator~.)+
  scale_fill_manual(values = color_mapping)
df_pp |> 
  group_by(predator, primary_producer, iso_week = isoweek(sample_week), station) |> 
  summarise(contrib = mean(rel_primary_production), .groups = "drop") |> 
  ggplot(aes(x = iso_week, y = contrib, col = predator, fill = predator))+
  geom_line(stat = "identity")+
  facet_grid(primary_producer~., scales = "free")

# PCoA + metadata helper
process_pcoa <- function(pcoa_result, metadata, id_cols = 1:3) {
  coords <- as.data.frame(pcoa_result$vectors[, 1:2])
  names(coords) <- c("Axis1", "Axis2")
  
  coords |>
    bind_cols(metadata[, id_cols]) |>
    mutate(
      iso_week = isoweek(sample_week),
      season = case_when(
        iso_week %in% 1:10   ~ "Winter",
        iso_week %in% 11:22  ~ "Spring",
        iso_week %in% 23:35  ~ "Summer",
        iso_week %in% 36:48  ~ "Autumn",
        iso_week %in% 49:53  ~ "Winter",
        TRUE ~ NA_character_
      ),
      month = month(sample_week),
      month_abb = factor(month.abb[month], levels = month.abb)
    )
}

# Prepare phytoplankton data (needs pivoting first)
plot_pcoa <- function(df, eig,envfit, title) {
  p1 <- ggplot(df, aes(x = Axis1, y = Axis2)) +
    
    geom_point(shape = 21, size = 1.5, mapping = aes(fill = predator)) +
    facet_wrap(~month_abb) +
    scale_fill_manual(values = color_mapping) +
    coord_fixed() +
    theme_bw() +

    labs(
      title = title,
      x = paste0("Axis 1 (", round(100 * eig[1], 1), "%)"),
      y = paste0("Axis 2 (", round(100 * eig[2], 1), "%)")
    )
  p2 <- ggplot(envfit, mapping = aes(x = Axis1, y = Axis2)) + 
    geom_segment(mapping = aes(x = Axis1, y = Axis2, yend = 0, xend = 0),  arrow = arrow(length = unit(2, 'mm'), ends = "first")) +
    scale_y_continuous(limits = c(-1,1))+
    scale_x_continuous(limits = c(-1,1))+
    ggrepel::geom_text_repel(mapping = aes(x = Axis1, y = Axis2, label = prey), size = 2) +
    coord_fixed() +
    theme_bw()+
    labs(
      x = paste0("Axis 1 (", round(100 * eig[1], 1), "%)"),
      y = paste0("Axis 2 (", round(100 * eig[2], 1), "%)")
    )
  p1 + p2 + plot_layout(guides = "collect", width = c(3,1), axis_titles = "collect")
}
# For phytoplankton
mat_pp <- df_pp |>
  pivot_wider(names_from = primary_producer, values_from = rel_primary_production, values_fill = 0)
bray_dist_pp <- vegdist(mat_pp[, 4:15], method = "bray")
set.seed(100)
pcoa_result_pp <- pcoa(bray_dist_pp)
site_scores_pp <- pcoa_result_pp$vectors
envfit_pp <- as.data.frame(scores(envfit(site_scores_pp, mat_pp[4:15], permutations = 999), display = c("vectors"))) |> 
  rownames_to_column(var = "prey") |> 
  rename("Axis1" = Axis.1, "Axis2" = Axis.2)
eig_pp <- pcoa_result_pp$values$Relative_eig
pcoa_df_pp <- process_pcoa(pcoa_result_pp, metadata = mat_zp)
pcoa_pp <- plot_pcoa(pcoa_df_pp, eig_pp,envfit_pp, "PCoA of Energy Fluxes based on relative contribution of phytoplankton to fish")

# For zooplankton
bray_dist_zp <- vegdist(mat_zp[, 4:10], method = "bray")
set.seed(100)
pcoa_result_zp <- pcoa(bray_dist_zp)
site_scores_zp <- pcoa_result_zp$vectors
envfit_zp <- as.data.frame(scores(envfit(site_scores_zp, mat_zp[4:10], permutations = 999), display = c("vectors"))) |> 
  rownames_to_column(var = "prey") |> 
  rename("Axis1" = Axis.1, "Axis2" = Axis.2)
eig_zp <- pcoa_result_zp$values$Relative_eig
pcoa_df_zp <- process_pcoa(pcoa_result_zp, metadata = mat_zp)
pcoa_zp <- plot_pcoa(pcoa_df_zp, eig_zp,envfit_zp, "PCoA of Energy Fluxes based on relative contribution of phytoplankton to fish")
pcoa_zp/pcoa_pp
ggsave("./output/figure/pcoa.pdf", dpi = 500, width = 10, height = 10)

# Function to generate a tidy table from pairwise Adonis results
ptable_pwadonis <- function(pwa) {
  pwa[-1] |> 
    imap_dfr(~ tibble(
      Pairs = .y,
      F.Model = .x$F[1],
      R2 = .x$R2[1],
      `Pr(>F)` = .x$`Pr(>F)`[1],
      Df = paste(.x$Df[1], .x$Df[2], sep = ";")
    )) |> 
    separate(Pairs, into = c("Group1", "Group2"), sep = "_vs_") |> 
    mutate(
      group_min = pmin(Group1, Group2),
      group_max = pmax(Group1, Group2)
    ) |> 
    select(-Group1, -Group2) |> 
    rename(Group1 = group_min, Group2 = group_max) |> 
    mutate(
      LabelA = str_c("Between ", Group1, " and ", Group2),
      LabelB = str_c("Between ", Group2, " and ", Group1)
    )
}

# Function to perform posthoc pairwise comparison for a given month
posthoc_pairwise_predator <- function(i, mat) {
  require(pairwiseAdonis)
  mat_month <- mat |>  filter(month == i)
  
  # Select all numeric columns except 'month' for distance matrix calculation
  dist_mat <- mat_month |> 
    select(where(is.numeric) & !matches("month")) |> 
    vegdist(method = "bray")
  
  # Perform pairwise adonis
  pwa <- pairwise.adonis2(dist_mat ~ predator, data = mat_month)
  ptable_pwadonis(pwa) |>  mutate(month = i)
}

# Function to run PERMANOVA and pairwise posthoc for a dataset
run_permanova_workflow <- function(mat, seed = 100) {
  set.seed(seed)
  mat <- mat |> 
    mutate(month = month(sample_week))
  
  # Global PERMANOVA using numeric columns except 'month'
  adonis_global <- adonis2(
    formula = vegdist(select(mat, where(is.numeric) & !matches("month")), method = "bray") ~ predator * month,
    data = mat,
    permutations = 999,
    strata = mat$month,
    parallel = 8
  )
  
  # Pairwise posthoc test per month
  pairwise_posthoc <- map_dfr(1:12, ~posthoc_pairwise_predator(.x, mat)) |> 
    group_by(month) |> 
    mutate(p_adj = p.adjust(`Pr(>F)`, method = "BH")) |> 
    ungroup()
  
  # Betadisper per month (with corrected filtering and group assignment)
  betadisper_results <- map_dfr(1:12, function(i) {
    disp_mat <- mat |> filter(month == i)
    
    dist_mat <- vegdist(select(disp_mat, where(is.numeric) & !matches("month")), method = "bray")
    disp <- betadisper(dist_mat, group = disp_mat$predator)
    disp_test <- permutest(disp, permutations = 999)
    
    tibble(
      month = i,
      F_value = disp_test$tab[1, "F"],
      P_value = disp_test$tab[1, "Pr(>F)"]
    )
  })
  
  list(
    global = adonis_global,
    pairwise = pairwise_posthoc,
    dispersion = betadisper_results
  )
}

# Example: Running the workflow for 'mat_zp'
mat_zp <- mat_zp |> 
  mutate(month = month(sample_week)) # Add 'month' column

# Running PERMANOVA and pairwise comparisons for 'mat_zp'
results_zp <- run_permanova_workflow(mat_zp)

# Example: Running the workflow for 'mat_pp'
mat_pp <- mat_pp |> 
  mutate(month = month(sample_week)) # Add 'month' column

# Running PERMANOVA and pairwise comparisons for 'mat_pp'
results_pp <- run_permanova_workflow(mat_pp)

results_pp$global
results_zp$global

results_pp$pairwise |> 
  mutate(type = "pp") |> 
  left_join(results_pp$dispersion |> rename("dispersionP" = P_value), by = join_by(month)) |> 
  rbind(results_zp$pairwise |> 
          mutate(type = "zp") |> 
  left_join(results_zp$dispersion |> rename("dispersionP" = P_value), by = join_by(month)) ) |> 
  left_join(tibble(month = 1:12, month_abb = month.abb), by = join_by(month)) |> 
  mutate(month_abb = factor(month_abb, levels = month.abb),
         sig = ifelse(p_adj <= 0.05, "*", ""),
         dispersion = ifelse(dispersionP <= 0.05, "!", "")) |> 
  
  ggplot(aes(x = Group1, y = Group2, fill = R2, label = paste(sig, dispersion))) +
  geom_tile(col = "black") +
  geom_text()+
  facet_grid(type~month_abb) +
  coord_fixed() +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  scale_fill_gradient2(high = "black", low = "#DDFDFE", mid = "#D16666", midpoint = 0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.grid = element_blank()) +
  labs(x = NULL, y = NULL)
ggsave("./output/figure/permanova.pdf", dpi = 500, height = 2.5, width = 10)
rm(bray_dist_pp, df_pp, fish_fluxes, mat_pp, mat_zp, pcoa_result_pp, rel_contribution_fluxes, zp_fluxes, tot_fluxes, site_scores_pp)

## Impact of forage ratios ----
# or how the forage ratios deviate from the null model
isoweek_null <- timeseries_null |> 
  group_by(predator, prey, iso_week = isoweek(sample_week), station) |> 
  summarise(avg_fluxes = mean(flux, na.rm = T), .groups = "drop")

# Reusable plot function
plot_flux_difference <- function(trophic_level_filter) {
  isoweek_fluxes |> 
    left_join(isoweek_null |> rename("null" = avg_fluxes),
              by = c("iso_week", "station", "predator", "prey")) |> 
    mutate(
      predator = factor(predator, levels = node_data$node_name),
      prey = factor(prey, levels = node_data$node_name),
      trophic_level = case_when(
        prey %in% node_data$node_name[node_data$trophic_level == 2] ~ "zooplankton",
        prey %in% node_data$node_name[node_data$trophic_level == 1] ~ "phytoplankton"
      )
    ) |>
    filter(trophic_level == trophic_level_filter,
           iso_week %in% 2:52) |>
    ggplot(aes(x = iso_week)) +
    
    geom_ribbon(aes(y = mean, ymin = lower, ymax = upper),
                alpha = .2, linewidth = .2,
                fill = "black", col = "black") +
    
    geom_line(aes(y = mean, color = "Selectivity")) +
    geom_line(aes(y = null, color = "No selectivity")) +
    
    facet_grid(
      prey ~ predator, 
      scales = "free_y"
    ) +
    
    scale_color_manual(
      name = NULL,
      values = c("Selectivity" = "black", "No selectivity" = "#ff7f00")
    ) +
    
    scale_x_continuous(
      breaks = seq(1, 52.1775, 4.348125*2),
      labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"), expand = c(0, 0)
    ) +
    
    theme_bw() +
    labs(
      y = "Fluxes \n [kJ/day/m2]",
      x = NULL
    )
}
phytoplankton_diff <- plot_flux_difference("phytoplankton")
zooplankton_diff   <- plot_flux_difference("zooplankton")
zooplankton_diff + phytoplankton_diff + plot_layout(width = c(3,8), guides = "collect", axis_titles = "collect")
## Over the entire timeseries ----
tot_null <- timeseries_null |> 
  group_by(station, predator, prey) |> 
  summarise(flux = mean(flux, na.rm = T), .groups = "drop")
flux_diff <- station_fluxes |> 
  left_join(tot_null, by = c("station", "predator", "prey")) |>
  mutate(
    sig = ifelse((lower - flux > 0) | (upper - flux < 0), "Significant", "Not significant"),
    prey = factor(prey, levels = node_data$node_name) |> fct_rev(),
    predator = factor(predator, levels = node_data$node_name),
    trophic_level = case_when(
      prey %in% node_data$node_name[node_data$trophic_level == 2] ~ "zooplankton",
      prey %in% node_data$node_name[node_data$trophic_level == 1] ~ "phytoplankton"
    )
  ) |> 
  group_by(station, predator, trophic_level, sig) |> 
  mutate(
    rel_mean  = (mean - flux) / flux * 100,
    rel_upper = (upper - flux) / flux * 100,
    rel_lower = (lower - flux) / flux * 100,
    rel_null  = 0,
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
    ggplot(aes(y = predator, x = rel_mean, xmax = rel_upper, xmin = rel_lower,
               fill = sig)) +
    geom_vline(xintercept = 0, color = "#388659", linewidth = 2) +
    geom_errorbarh(height = 0.2, position = position_dodge2(width = 0.2, padding = 3)) +
    geom_point(shape = 21, size = 2, position = position_dodge2(width = 0.2, padding = 3)) +
    facet_grid(. ~ prey, scales = "free") +
    scale_fill_manual(values = c("black", "white")) +
    labs(x = "Difference from null model [%]", y = NULL) +
    theme_bw()
}

plot_zoop <- plot_flux_diff(flux_diff, "zooplankton")
plot_phyto <- plot_flux_diff(flux_diff, "phytoplankton")
(plot_zoop + ggtitle("Percentage difference with selectivity compared to presence absence over the entire timeseries")) / plot_phyto + plot_layout(heights = c(3,7), guides = "collect") 


# Predation pressure ----
# Join and calculate predation pressure
timeseries_pressure <- timeseries_nodes |> 
  filter(
    !name %in% c("Clupea", "Gasterosteus", "Sprattus")
  ) |> 
  rename("prey" = name) |> 
  select(prey, sample_week, biomass) |> 
  right_join(timeseries_fluxes, by = c("prey", "sample_week"), relationship = "many-to-many") |> 
  mutate(
    iso_week = isoweek(sample_week)
  ) |>
  group_by(prey, predator, sample_week, station, iso_week, year, biomass) |> 
  summarise(pressure_sum = sum(mean, na.rm = TRUE), .groups = "drop_last") |> 
  summarise(pressure = sum(pressure_sum, na.rm = TRUE) / sum(biomass, na.rm = TRUE), .groups = "drop") |> 
  filter(iso_week %in% 2:52) |> 
  mutate(
    trophic_level = case_when(
      prey %in% node_data$node_name[node_data$trophic_level == 2] ~ "zooplankton",
      prey %in% node_data$node_name[node_data$trophic_level == 1] ~ "phytoplankton"
    )
  ) |> 
  pivot_wider(names_from = prey, values_from = pressure, values_fill = 0) |>
  pivot_longer(7:25, names_to = "prey", values_to = "pressure") |> 
  group_by(predator, prey) |> 
  filter(sum(pressure)>0) |> 
  ungroup()

timeseries_pressure |> 
  group_by(iso_week, predator, prey) |> 
  summarise(pressure = mean(pressure), .groups = "drop") |> 
  mutate(prey = factor(prey, levels = node_data$node_name),
         predator = factor(predator, levels = node_data$node_name)) |> 
  ggplot(aes(x = iso_week, y = pressure, fill =  prey))+
  geom_area()+
  facet_grid(predator~., scales = "free") +
  scale_fill_manual(values = c("Total" = "black",color_mapping))+
  labs(y = "Incoming fluxes / prey biomass [kJ/day/g]",
       x = NULL)+    
  scale_x_continuous(breaks = seq(1, 52.1775, 4.348125),
                     labels = month.abb, expand = c(0, 0)) +
  theme_bw()

timeseries_pressure |> 
  group_by(iso_week, predator, prey) |> 
  summarise(pressure = mean(pressure), .groups = "drop") |> 
  mutate(prey = factor(prey, levels = node_data$node_name),
         predator = factor(predator, levels = node_data$node_name)) |> 
  ggplot(aes(x = iso_week, y = pressure, fill =  predator))+
  geom_area()+
  facet_wrap(.~prey, scales = "free", ncol = 3) +
  scale_fill_manual(values = c("Total" = "black",color_mapping))+
  labs(y = "Outgoing fluxes / biomass [kJ/day/g]",
       x = NULL)+    
  scale_x_continuous(breaks = seq(1, 52.1775, 4.348125),
                     labels = month.abb, expand = c(0, 0)) +
  theme_bw()
  
# Work in progress ----
# Join and calculate predation pressure
weekly_biomasses <- read_csv("./data/processed/interpolation/weekly_biomasses.csv")
timeseries_pred_pressure <- weekly_biomasses |> 
  filter(
    station_name == "BY31 LANDSORTSDJ",
    !node_name %in% c("Clupea", "Gasterosteus", "Sprattus")
  ) |> 
  rename_with(~c("prey", "station"), .cols = c("node_name", "station_name")) |> 
  select(prey, station, sample_week, biomass) |> 
  right_join(timeseries_fluxes, by = c("prey", "station", "sample_week"), relationship = "many-to-many") |> 
  mutate(
    iso_week = isoweek(sample_week),
    year = as.factor(year(sample_week))
  ) |>
  group_by(prey, predator, sample_week, station, iso_week, year) |> 
  summarise(pressure = sum(mean, na.rm = TRUE) / sum(biomass, na.rm = TRUE), .groups = "drop") |> 
  filter(iso_week %in% 2:52) |> 
  mutate(
    pressure = replace_na(pressure, 0),
    trophic_level = case_when(
      prey %in% node_data$node_name[node_data$trophic_level == 2] ~ "zooplankton",
      prey %in% node_data$node_name[node_data$trophic_level == 1] ~ "phytoplankton"
    )
  ) |> 
  right_join(
weekly_biomasses |> 
  filter(station_name == "BY31 LANDSORTSDJ") |> 
  rename("predator" = node_name,
         "station" = station_name) |> 
  select(sample_week, predator, station, biomass),
by = c("sample_week", "predator", "station")) |> 
  filter(biomass > 0) |> 
  mutate(norm_pressure = pressure/biomass,
         norm_pressure = replace_na(norm_pressure, 0))
# Add total predation pressure across prey
timeseries_pred_pressure_with_tot <- timeseries_pred_pressure |> 
  group_by(prey = "Total",predator, sample_week, station, iso_week, year, trophic_level) |> 
  summarise(pressure = sum(pressure, na.rm = TRUE),
            norm_pressure = sum(norm_pressure, na.rm = TRUE),.groups = "drop") |> 
  bind_rows(timeseries_pred_pressure) |> 
  mutate(prey = factor(prey, levels = c("Total", node_data$node_name)),
         predator = factor(predator, levels = c(node_data$node_name)),
         norm_pressure = replace_na(norm_pressure, 0)) |> filter(prey == "Total",
                                                                                 !is.na(prey),
                                                                                 !is.na(norm_pressure),
                                                                                 !is.na(trophic_level))

# Aggregate by isoweek
isoweek_pred_pressure <- timeseries_pred_pressure_with_tot |> 
  group_by(prey, predator, iso_week, station, trophic_level) |> 
  summarise(median_pressure = median(pressure, na.rm = TRUE),
            pressure_upper = quantile(pressure, 0.75, na.rm = TRUE),
            pressure_lower = quantile(pressure, 0.25, na.rm = TRUE),
            median_norm_pressure = median(norm_pressure, na.rm = TRUE),
            pressure_norm_upper = quantile(norm_pressure, 0.75, na.rm = TRUE),
            pressure_norm_lower = quantile(norm_pressure, 0.25, na.rm = TRUE),
            .groups = "drop") |> 
  na.omit()

ggplot() +
    geom_point(data = timeseries_pred_pressure_with_tot,
            mapping = aes(x = iso_week, y = norm_pressure, group = year), col = "black", size = .2)+
  geom_ribbon(alpha = .1,
              data = isoweek_pred_pressure,mapping = aes(x = iso_week, y = median_norm_pressure, ymax = pressure_norm_upper, ymin = pressure_norm_lower,  col = predator, fill = predator))+
  geom_line(linewidth = 1,
            data = isoweek_pred_pressure,mapping = aes(x = iso_week, y = median_norm_pressure, col = predator)) +
  
  facet_wrap(predator+trophic_level~. , scales = "free_y", nrow = 3) +
  theme_bw() +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb, expand = c(0, 0)
  ) +
  scale_color_manual(values = color_mapping)+
  scale_fill_manual(values = color_mapping) +
  #scale_y_log10()+
     labs(
    x = NULL,
    y = "Predation pressure \n (median +- interquartile ranges) \n [kJ/day/g / g of predator]"
  )

ggplot(data = timeseries_pred_pressure_with_tot,
       mapping = aes(x = sample_week, y = norm_pressure, group = year)) +
  geom_line(col = "black", size = .2)+
  geom_smooth(
    method = "loess", se = FALSE, color = "red", linewidth = 0.8,
    aes(group = 1)  # smooth across all years per facet
  )+
  facet_wrap(predator+trophic_level~. , scales = "free_y", nrow = 3) +
  theme_bw() +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb, expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Predation pressure \n (median +- interquartile ranges) \n [kJ/day/g / g of predator]"
  )




timeseries_fluxes |> 
  left_join(node_data, by = c(predator = "node_name")) |> 
  group_by(sample_week, iso_week, trophic_level, station) |> 
  summarise(flux = sum(mean, na.rm = T), .groups = "drop") |>
  mutate(interaction = ifelse(trophic_level == 2, "from phytoplankton to zooplankton", "from zooplankton to fish"),
         interaction = factor(interaction, levels = c("from zooplankton to fish","from phytoplankton to zooplankton"))) |>
  filter(year(sample_week)>2007) |> 
  group_by(interaction,year = year(sample_week), station) |> 
  summarise(flux = mean(flux, na.rm = T)* 365, .groups = "drop") |> 
  ggplot(aes(x=year, y = flux, col = interaction))+
  geom_line() +
  geom_point()+
  labs(y = "annual fluxes \n [kJ/yr/m2]")
timeseries_fluxes |> 
  left_join(node_data, by = c(predator = "node_name")) |> 
  group_by(sample_week, trophic_level, station, type = "selectivity") |> 
  summarise(flux = sum(mean, na.rm = T), .groups = "drop") |>
  mutate(interaction = ifelse(trophic_level == 2, "from phytoplankton to zooplankton", "from zooplankton to fish"),
         interaction = factor(interaction, levels = c("from zooplankton to fish","from phytoplankton to zooplankton"))) |> 
  ggplot(aes(x = sample_week, y = flux))+
  geom_line()+
  facet_grid(interaction~., scale = "free") +
  theme_bw()
timeseries_fluxes |> 
  left_join(node_data, by = c(predator = "node_name")) |> 
  group_by(sample_week, trophic_level, station, type = "selectivity") |> 
  summarise(flux = sum(mean, na.rm = T), .groups = "drop") |>
  mutate(interaction = ifelse(trophic_level == 2, "pp_zp", "zp_fish")) |> 
  select(-trophic_level) |>
  pivot_wider(names_from = interaction, values_from = flux) |> 
  mutate(efficiency = zp_fish/pp_zp) |> 
  ggplot(aes(x = sample_week, y = efficiency, col = factor(month(sample_week))))+
  geom_line(col = "black")+
#  geom_point()+
  theme_bw()+
  labs(y = "fluxes ratio \n zp to fish / pp to zp ")

