#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))


#Prepare data ----
# Check if biomass data exists, otherwise create it.
if(!file.exists(file.path("data", "processed", "interpolation", "weekly_biomasses.csv"))) {
  system(paste("nohup Rscript", file.path("code", "InterpolateWeekly.R")))
}

# Check if predator selectivity data exists, otherwise create it.
if(!file.exists(file.path("data", "processed", "predator_selectivity.csv"))) {
  system(paste("nohup Rscript", file.path("code", "CombineMetabarcoding.R")))
}

# Merge selectivity and biomass ----
biomass <-
  read_csv(file.path("data", "processed", "interpolation", "weekly_biomasses.csv")) |> 
  group_by(sample_week, station_name) |>
  mutate(rel_biomass = biomass / sum(biomass)) |> 
  select(node_prey = node_name, sample_week, station_name, biomass)
ForageRatios <- 
  read_csv(file.path("data", "processed", "predator_selectivity.csv")) |> 
  ungroup() |>
  mutate(rra_gut = replace_na(rra_gut, 0),
         rra_env = replace_na(rra_env, 0),
         ForageRatio = rra_gut/rra_env) |>
  filter(rra_env > 0) |> 
  left_join(biomass) |> 
  group_by(sample_id) |> 
  mutate(rel_biomass = biomass / sum(biomass)) |> 
  ungroup()
# Quick visualisation
ggplot(ForageRatios |> filter(rel_biomass>0), aes(x = rel_biomass, y = ForageRatio+1))+
  geom_point()+
  scale_y_log10()+
  facet_grid(node_prey~node_predator)
# Calculate the average forage ratio for each predator-prey pairs
average_forage_ratios <-
  ForageRatios |> 
  group_by(node_predator, node_prey) |> 
  summarise(ForageRatio = mean(ForageRatio, na.rm = T), .groups = "drop")

# Estimate the forage ratio based on a type II response:
# Step 1: Create a bootstrap function
bootstrap_fit <- function(df, model, coef, prey, predator) {
  df_boot <- df[sample(nrow(df), replace = TRUE), ]  # Resample with replacement
  
  boot_fit <- tryCatch(
    nlsLM(model,
          start = list(a = coef[1, 1], h = coef[2, 1]),
          data = df_boot,
          lower = c(a_min, h_min),
          upper = c(a_max, h_max),
          control = nls.lm.control(maxiter = 500)),
    error = function(e) NULL  
  )
  
  if (!is.null(boot_fit)) {
    boot_coef <- summary(boot_fit)$coefficients
    return(tibble(
      node_predator = predator,
      node_prey = prey,
      a = boot_coef[1, 1],
      h = boot_coef[2, 1]
    ))
  } else {
    return(NULL) 
  }
}

# As the previous function can return NULL values, we want to make sure that we have 1000 values that are non-NULL 
run_bootstraps <- function(df, model, coef, prey, predator, n_boot = 1000) {
  successful_results <- list()
  attempts <- 0
  
  while (length(successful_results) < n_boot) {
    result <- bootstrap_fit(df, model, coef, prey, predator)
    attempts <- attempts + 1
    
    if (!is.null(result)) {
      successful_results[[length(successful_results) + 1]] <- result
    }
    
    # Optional: break infinite loop (e.g., max 10x more attempts than desired)
    if (attempts > n_boot * 10) {
      warning(glue::glue("Reached max attempts for {predator}-{prey}: {length(successful_results)} fits collected."))
      break
    }
  }
  
  bind_rows(successful_results) %>%
    mutate(Iteration = row_number())
}

# Step 2: set up everything that will be needed
# Libraries and start a multisession so it paralellise
suppressPackageStartupMessages(library(furrr))
plan(multisession) 
library(minpack.lm)
# parameters that are used in the model
set.seed(100)
n_boot = 1200
a_values = c(0, 1, 10, 100, 1000) ; a_min = 0 ; a_max = 1000
h_values = c(1, 60, 60*24) ; h_min = 0 ; h_max = 60*60*24
# the model itself
model <- ForageRatio ~ a  / (1 + a * h * rel_biomass)
# empty tibble and list for storing output
results <- tibble()  # Stores best-fit parameters
boot_results <- list()  # Stores bootstrapping results

# Step 3: run the model through a for loop
# Split data by Barcode, and Species_English
df_list <-
  ForageRatios |>
  group_split(node_predator, node_prey, .keep = TRUE)

# Fit models and perform bootstrapping in parallel
for (df in df_list) {
  if (nrow(df) <= 5) next  # Skip if there are 5 or fewer observations
  
  prey <- unique(df$node_prey)
  predator <- unique(df$node_predator)
  
  best_fit <- NULL
  best_residual_dispersion <- Inf  
  
  for (a_start in a_values) {
    for (h_start in h_values) {
      
      fit_attempt <- tryCatch(
        nlsLM(model,
              start = list(a = a_start, h = h_start),
              data = df,
              lower = c(a_min, h_min),
              upper = c(a_max, h_max),
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
    plot(fitted(best_fit), resid(best_fit), main = paste(predator, prey, sep = "-"))
    abline(h = 0)
    # Now store
    coef <- summary(best_fit)$coefficients
    
    results <-
      results |>
      bind_rows(
        tibble(node_predator = predator,
               node_prey = prey,
               a = coef[1, 1],
               h = coef[2, 1]
        )
      )
    
    # Run bootstrapping in parallel
#    boot_results[[paste(prey, predator)]] <- future_map_dfr(
 #     1:n_boot, 
  #    ~ bootstrap_fit(df, model, coef, prey, predator), 
   #   .options = furrr_options(seed = TRUE)# For reproducibility
    #) |> 
     # mutate(Iteration = row_number()) #|> 
      #filter(Iteration %in% 1:1000)
    boot_results[[paste(prey, predator)]] <- run_bootstraps(df, model, coef, prey, predator, n_boot = 1000)
  }
}

# Summarise the confidence intervals
bootstrapped_values <-
  bind_rows(boot_results) |> 
  group_by(node_prey, node_predator) |> 
  summarise(
    AVG_a = median(a),
    AVG_h = median(h),
    LOW_a = quantile(a, 0.025),
    LOW_h = quantile(h, 0.025),
    HIGH_a = quantile(a, 0.975),
    HIGH_h = quantile(h, 0.975),
    .groups = "drop"
  ) |> 
  mutate(
    a = paste0(round(AVG_a, 2), " [", round(LOW_a, 2), ";", round(HIGH_a, 2), "]"),
    h = paste0(round(AVG_h, 2), " [", round(LOW_h, 2), ";", round(HIGH_h, 2), "]")
  ) |> 
  select( node_predator,node_prey, a, h)|> 
  arrange(node_predator, node_prey)
print(bootstrapped_values, n = nrow(bootstrapped_values))
# Step 4: Calculate confidence interval:
split_boot <- bind_rows(boot_results) |> group_split(Iteration, node_predator)
read_csv(file.path("data", "processed", "interpolation", "weekly_biomasses.csv")) |>
  filter(station_name == "BY31 LANDSORTSDJ") |> 
  mutate(type = ifelse(node_name %in% c("Sprattus", "Clupea", "Gasterosteus"), "fish", ifelse(node_name %in% c("Acartia", "Bosmina", "Centropages", "Eurytemora", "Evadne", "Pseudocalanus", "Synchaeta", "Temora"), "zooplankton", "phytoplankton")))
  group_by(sample_week, station_name) |>
  mutate(rel_biomass = biomass / sum(biomass)) |> 
  ggplot(aes(x= sample_week, y = rel_biomass, fill = node_name))+
  geom_area(stat = "identity")
###############################################################################
modelled_forage_ratios <-
  bind_rows(boot_results) |> 
  group_by(node_prey, node_predator) |> 
  summarise(
    a = median(a),
    h = median(h),

    .groups = "drop"
  ) |> 
  right_join(average_forage_ratios) |> 
  arrange(node_predator, node_prey)
rel_biomass_seq <- seq(0.001, 1, length.out = 100)
modelled_forage_ratios |> 
  cross_join(tibble(rel_biomass = rel_biomass_seq)) |> 
  mutate(projected_ForageRatio = ifelse(!is.na(a) & !is.na(h), ((a) / (1 + a * h * rel_biomass)), ForageRatio)) |> 
  ggplot(aes(x = rel_biomass, y = projected_ForageRatio))+
  geom_line()+
  facet_grid(node_predator~node_prey, scales = "free_y")


