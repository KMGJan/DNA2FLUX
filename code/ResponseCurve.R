#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(minpack.lm))


# this will make sure that the interpolated zooplankton data exists: 
if(!file.exists(file.path("data", "processed", "interpolation", "zooplankton_biomass.csv"))) {
  getinterpolation <- file.path("code", "InterpolateWeekly.R")
  # Get the monitoring data
  system(paste("nohup Rscript", getinterpolation))
  rm(getinterpolation)
} else {
  cat(paste(
    file.path("data", "processed", "interpolation", "zooplankton_biomass.csv"),
    "exists. Proceed with the analyses...")
    )
}

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
  #  Keep genus that are in the zooplankton dataset
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
  # Keep genus that are in the zooplankton dataset
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
# Quick inspection of the data ----------------------------------------------- #
df <-
  df_COI |> 
  bind_rows(df_18S) |> 
  select(-c("Supergroup", "Division", "Subdivision"))

df |> 
  group_by(sample_ID, Taxa, barcode, organism) |>
  summarise(RRA = sum(RRA, na.rm = T), .groups = "drop") |> 
  ggplot(aes(x = sample_ID, y = RRA, fill = Taxa)) +
  geom_bar(stat = "identity") +
  facet_grid(barcode ~ organism, scales = "free")
# Combine fish and environment ----------------------------------------------- #
indices_df <-
  # Combine the fish samples
  filter(df, organism != "WP2") |> 
  # With the environmental samples
  left_join(filter(df, organism == "WP2"),
            by = c("collection_date", "env_broad_scale", "env_medium", "station_ID", "barcode", "ASV_ID", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "sequence", "Taxa"),
            suffix = c(".gut", ".env"),
            relationship = "many-to-many" # Because 1 environmental sample matches up to 9 gut samples
  ) |> 
  #Remove sequences that are not present in the water
  distinct() |> 
  # Rename some key parameters for better readability
  rename("sample_ID" = sample_ID.gut,
         "organism" = organism.gut,
         "lat_lon" = lat_lon.gut) |>
  # Rgut and Renv correspond to all reads from the same genus
  group_by(sample_ID, station_ID, collection_date, organism, barcode, Taxa, lat_lon) |> 
  summarise(Rgut = sum(RRA.gut, na.rm = T),
            Renv = sum(RRA.env, na.rm = T),
            .groups = "drop") |>
  # Remove Renv < 0 as it is impossible to calculate a ratio for these values
  filter(Renv > 0, Taxa != "Other") |> 
  # Calculate the forage ratios
  group_by(sample_ID, barcode) |> 
  mutate(Rgut = Rgut / sum(Rgut, na.rm = T),
         Renv = Renv / sum(Renv, na.rm = T),
         ForageRatio = Rgut / Renv,
         sample_week = floor_date(collection_date, unit = "week", week_start = 1)) |>    # Change date to the first day of the week so it will match the zooplankton dataset
  ungroup()
# Quick inspection of the data ----------------------------------------------- #
indices_df |>
  ggplot(aes(x = organism, y = ForageRatio+1, fill = Taxa))+
  geom_boxplot() +
  labs(y = "Forage ratio (Rgut/Renv)")+
  scale_y_log10() +
  facet_grid(barcode~organism, scales = "free")+
  geom_hline(yintercept = 2)+
  theme_bw()+
  theme(axis.text.x = element_blank())
# Combine with monitoring dataset ----------------------------------------------
test_df <-
  read_csv(file.path("data", "processed", "interpolation", "zooplankton_biomass.csv"), show_col_types = F) |> 
  #work in relative biomass
  group_by(sample_week) |> 
  mutate(rel_biomass = (biomass / sum(biomass, na.rm = T)),
         Taxa = node_name) |> 
  ungroup() |> 
  right_join(indices_df, 
             by = c("sample_week", "Taxa"),
             relationship = "many-to-many") |> 
  filter(Taxa != "Other")
test_df |> 
  filter(barcode != "Both") |> 
  ggplot(aes(x = rel_biomass, y = ForageRatio))+
  geom_point(mapping = aes(col = barcode))+
  facet_wrap(Taxa~organism, scales = "free", ncol = 3)+
  geom_hline(yintercept = 1)
# Attempt to fit the model with confidence intervals----------------------------
# Set up parallel processing
suppressPackageStartupMessages(library(furrr))
plan(multisession) 
# The model
model <- ForageRatio ~ (a * rel_biomass^k) / (1 + a * h * rel_biomass^k) / rel_biomass

# Function for bootstrapping
bootstrap_fit <- function(df, model, coef, genus, fish) {
  df_boot <- df[sample(nrow(df), replace = TRUE), ]  # Resample with replacement
  
  boot_fit <- tryCatch(
    nlsLM(model,
          start = list(a = coef[1, 1], h = coef[2, 1], k = coef[3, 1]),
          data = df_boot,
          lower = c(a_min, h_min, k_min),
          upper = c(a_max, h_max, k_max),
          control = nls.lm.control(maxiter = 100)),
    error = function(e) NULL  
  )
  
  if (!is.null(boot_fit)) {
    boot_coef <- summary(boot_fit)$coefficients
    return(tibble(
      organism = fish,
      Taxa = genus,
      a = boot_coef[1, 1],
      h = boot_coef[2, 1],
      k = boot_coef[3, 1]
    ))
  } else {
    return(NULL) 
  }
}

# Some parameters that are used later on:
set.seed(100)
n_boot = 1200
a_values = c(0, 1, 10, 100, 1000) ; a_min = 0 ; a_max = 1000
h_values = c(1, 60, 60*24) ; h_min = 1 ; h_max = 60*60*24
k_start = 1 ; k_min = 1 ; k_max = 1
# Run the model to estimate a, h and k -----------------------------------------
# fit_me <- function(your_df){  

# Split data by Barcode, Genus, and Species_English
df_list <-
  test_df |>
  group_split(Taxa, organism, .keep = TRUE)
# Initialize results storage
results <- tibble()  # Stores best-fit parameters
boot_results <- list()  # Stores bootstrapping results

# Fit models and perform bootstrapping in parallel
for (df in df_list) {
  genus <- unique(df$Taxa)
  fish <- unique(df$organism)
  
  best_fit <- NULL
  best_residual_dispersion <- Inf  
  
  for (a_start in a_values) {
    for (h_start in h_values) {
      
      fit_attempt <- tryCatch(
        nlsLM(model,
              start = list(a = a_start, h = h_start, k = k_start),
              data = df,
              lower = c(a_min, h_min, k_min),
              upper = c(a_max, h_max, k_max),
              control = nls.lm.control(maxiter = 500)),
        error = function(e) NULL)
      
      
      if (!is.null(fit_attempt)) {
        residual_dispersion <- sum(residuals(fit_attempt)^2)
        
        if (!is.null(fit_attempt) && residual_dispersion < best_residual_dispersion) {
          best_residual_dispersion <- residual_dispersion
          best_fit <- fit_attempt
        }
        
      }
    }
  }
  
  # Store best-fit model results (without bootstrapping)
  if (!is.null(best_fit)) {
    # Residual visualisation
    plot(fitted(best_fit), resid(best_fit), main = paste(genus, fish, sep = "-"))
    abline(h = 0)
    # Now store
    coef <- summary(best_fit)$coefficients
    
    results <-
      results |>
      bind_rows(
        tibble(organism = fish,
               Taxa = genus,
               a = coef[1, 1],
               h = coef[2, 1],
               k = coef[3, 1]
        )
      )
    
    # Run bootstrapping in parallel
    boot_results[[paste(genus, fish)]] <- future_map_dfr(
      1:n_boot, 
      ~ bootstrap_fit(df, model, coef, genus, fish), 
      .options = furrr_options(seed = TRUE)
    ) |> 
      mutate(Iteration = row_number()) |> 
      filter(Iteration %in% 1:1000)
  }
}
# Summarise the confidence intervals
bootstrapped_values <-
  bind_rows(boot_results) |> 
  group_by(organism, Taxa) |> 
  summarise(
    AVG_a = mean(a),
    AVG_h = mean(h),
    AVG_k = mean(k),
    LOW_a = quantile(a, 0.025),
    LOW_h = quantile(h, 0.025),
    HIGH_a = quantile(a, 0.975),
    HIGH_h = quantile(h, 0.975),
    LOW_k = quantile(k, 0.025),
    HIGH_k = quantile(k, 0.975),
    .groups = "drop"
  ) |> 
  mutate(
    a = paste0(round(AVG_a, 2), " [", round(LOW_a, 2), ";", round(HIGH_a, 2), "]"),
    h = paste0(round(AVG_h, 2), " [", round(LOW_h, 2), ";", round(HIGH_h, 2), "]"),
    k = paste0(round(AVG_k, 2), " [", round(LOW_k, 2), ";", round(HIGH_k, 2), "]")
  ) |> 
  select(organism, Taxa, a, h, k)
print(bootstrapped_values, n = nrow(bootstrapped_values))
#
split_boot <- bind_rows(boot_results) |> group_split(Iteration, organism)
weekly_plankton <- read_csv(file.path("data", "processed", "interpolation", "zooplankton_biomass.csv"), show_col_types = F) |> 
  #work in relative biomass
  group_by(sample_week) |> 
  mutate(rel_biomass = (biomass / sum(biomass, na.rm = T)),
         Taxa = node_name) |> ungroup()
projectionForageRatio <-
  future_map(1:length(split_boot), ~{
    boot_data <- split_boot[[.x]]  # Extract each group of bootstrapped data
  
    weekly_plankton |> 
      left_join(boot_data, by = "Taxa", relationship = "many-to-many") |> 
      na.omit() |> 
      mutate(ForageRatio = (a * rel_biomass^k) / (1 + a * h * rel_biomass^k) / rel_biomass) |> 
      group_by(sample_week) |> 
      mutate(Prop = ForageRatio * biomass / sum(ForageRatio * biomass, na.rm = T),
             Chesson = ForageRatio / sum(ForageRatio, na.rm = T)) |> 
      ungroup()
    }, .options = furrr_options(seed = TRUE))  # Ensures reproducibility

projectionUncertainty <-
  bind_rows(projectionForageRatio) |>
  group_by(sample_week, organism, Taxa) |> 
  summarise(Med_prop = median(Prop, na.rm = T),
            Low_prop = quantile(Prop, 0.025, na.rm = T),
            High_prop = quantile(Prop, 0.975, na.rm = T),
            Med_fr = median(ForageRatio, na.rm = T),
            Low_fr = quantile(ForageRatio, 0.025, na.rm = T),
            High_fr = quantile(ForageRatio, 0.975, na.rm = T),
            Med_ch = median(Chesson, na.rm = T),
            Low_ch = quantile(Chesson, 0.025, na.rm = T),
            High_ch = quantile(Chesson, 0.975, na.rm = T),
            .groups = "drop")

rm(a_max, a_min, a_values, h_max, h_min, h_values, indices_df, k_max, k_min, k_start, n_boot, h_start, a_start, best_fit, best_residual_dispersion, coef, fish, genus)
rel_biomass_seq <- seq(0, 1, length.out = 100)

boot_prediction <-
  bind_rows(boot_results) |> 
  cross_join(tibble(rel_biomass = rel_biomass_seq)) |> 
  mutate(ForageRatio = (a * rel_biomass^k) / (1 + a * h * rel_biomass^k) / rel_biomass,
         Bgut = ForageRatio * rel_biomass) |> 
  group_by(Taxa, organism, rel_biomass) |> 
  summarise(fr_lower = quantile(ForageRatio, 0.025, na.rm = TRUE),
            fr_upper = quantile(ForageRatio, 0.975, na.rm = TRUE),
            Bgut_lower = quantile(Bgut, 0.025, na.rm = TRUE),
            Bgut_upper = quantile(Bgut, 0.975, na.rm = TRUE),
            .groups = "drop") 
# Quick inspection of the data ----------------------------------------------- #
results |> 
  na.omit() |> 
  cross_join(tibble(rel_biomass = rel_biomass_seq)) |> 
  mutate(ForageRatio = (a * rel_biomass^k) / (1 + a * h * rel_biomass^k) / rel_biomass,
         Bgut = ForageRatio * rel_biomass) |> 
  ggplot()+
  geom_ribbon(data = boot_prediction,
              mapping = aes(ymin = fr_lower, ymax = fr_upper, x = rel_biomass), alpha = .8, col = "black", linetype = 2)+
  
  geom_point(data = test_df |> filter(barcode != "Both"),
             mapping = aes(x = rel_biomass, y = ForageRatio), color = "black", shape = 21, alpha = .5, size = 1)+
  geom_line(linewidth = .5,
            mapping = aes(x = rel_biomass, y = ForageRatio))+
  
  facet_wrap(Taxa~organism, scales = "free", ncol = 3) +
  
  geom_hline(yintercept = 1)+
  theme_bw()+
  scale_y_log10()+
  labs(x = "Relative Biomass", y = "Forage ratio")

results |> 
  na.omit() |> 
  cross_join(tibble(rel_biomass = rel_biomass_seq)) |> 
  mutate(ForageRatio = (a * rel_biomass^k) / (1 + a * h * rel_biomass^k) / rel_biomass,
         Bgut = ForageRatio * rel_biomass) |> 
  ggplot()+
  geom_ribbon(data = boot_prediction,
              mapping = aes(ymin = fr_lower, ymax = fr_upper, x = rel_biomass), alpha = .8, col = "black", linetype = 2)+
  
  #geom_point(data = test_df |> filter(barcode != "Both"),
   #          mapping = aes(x = rel_biomass, y = ForageRatio), color = "black", shape = 21, alpha = .5, size = 1)+
  geom_line(linewidth = .5,
            mapping = aes(x = rel_biomass, y = ForageRatio))+
  
  facet_grid(Taxa~organism, scales = "free_y") +
  
  geom_hline(yintercept = 1)+
  theme_bw()+
  labs(x = "Relative Biomass", y = "Forage ratio")
results |> 
  na.omit() |> 
  cross_join(tibble(rel_biomass = rel_biomass_seq)) |> 
  mutate(ForageRatio = (a * rel_biomass^k) / (1 + a * h * rel_biomass^k) / rel_biomass,
         Bgut = ForageRatio * rel_biomass) |> 
  ggplot()+
  geom_ribbon(data = boot_prediction,
              mapping = aes(ymin = Bgut_lower, ymax = Bgut_upper, x = rel_biomass), alpha = .5, linetype = 2)+
  
  geom_line(mapping = aes(x = rel_biomass, y = Bgut),linewidth = .5) +
  
  
  facet_grid(organism~Taxa)+
  theme_bw()+
  coord_fixed()+
  labs(x = "Relative Biomass", y = "Predicted gut content")
# ---------------------------------------------------------------------------- #
# Projections
# ---------------------------------------------------------------------------- #
modeled_out <-
  weekly_plankton |> 
  left_join(results, by = "Taxa", relationship = "many-to-many") |>
  mutate(ForageRatio = (a*rel_biomass^k)/(1+a*h*rel_biomass^k)/rel_biomass) |> 
  group_by(sample_week, organism) |> 
  mutate(Prop = ForageRatio * biomass / sum(ForageRatio * biomass, na.rm = T)) |> ungroup() 
ggplot(mapping = aes(x = sample_week, col = organism, fill = organism))+
  geom_ribbon(data = projectionUncertainty,
              mapping = aes(ymin = Low_fr, ymax = High_fr), 
              alpha = .5, linewidth = .3)+
  geom_line(data = projectionUncertainty,
            mapping = aes(y = Med_fr),
            linewidth = .5)+
  facet_grid(Taxa~., scales = "free") +
  #scale_y_log10()+
  theme_bw()+
  geom_hline(yintercept = 1)+
  scale_color_manual(values = c(2,3,4))+
  labs(y = "Predicted Forage ratio")  
ggplot(mapping = aes(x = sample_week, col = Taxa, fill = Taxa))+
  geom_ribbon(data = projectionUncertainty,
              mapping = aes(ymin = Low_prop, ymax = High_prop), 
              alpha = .5, linewidth = .3)+
  geom_line(data = projectionUncertainty,
            mapping = aes(y = Med_prop),
            linewidth = .5)+
  facet_grid(organism~., scales = "free") +
  theme_bw()+
  labs(y = "Predicted link weight Proportion (Wij)")
projectionUncertainty |> 
  group_by(month(sample_week), organism, Taxa) |> 
  summarise(AVG = mean(Med_prop, na.rm = T)) |> ungroup() |> 
  na.omit() |> 
  ggplot(aes(x = as.factor(`month(sample_week)`), y = AVG, fill = organism)) +
  geom_bar(stat = "identity", alpha = 1, position = "dodge") +
  geom_errorbar(data = projectionUncertainty |> 
                  group_by(month(sample_week), organism, Taxa) |>
                  summarise(AVG_min = mean(Low_prop, na.rm = T),
                            AVG = mean(Med_prop, na.rm = T),
                            AVG_max = mean(High_prop, na.rm = T)) |> ungroup(),
                mapping = aes(ymin = AVG_min, ymax = AVG_max),position = "dodge") +
  facet_grid(Taxa~.) +
  theme_bw() +

  labs(y = "Predicted link weight Proportion (Wij)")
modeled_out |>
  mutate(Month = month(sample_week)) |> 
  na.omit() |> 
  group_by(Month, organism) |> 
  mutate(Prop = ForageRatio * biomass / sum(ForageRatio * biomass, na.rm =T)) |> ungroup() |> 
  
  ggplot(aes(x = as.factor(Month), y = Prop, fill = Taxa)) +
  geom_bar(stat = "identity", alpha = 1) +
  #  geom_line(aes(y = Biomass_g.m2/5, col = Genus))+
  facet_grid(organism~., scales = "free") +

  theme_bw() +
  labs(y = "Biomass contribution (Wij)")

