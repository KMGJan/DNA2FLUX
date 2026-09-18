#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(minpack.lm)
  library(furrr)
  library(progressr)
  library(patchwork)
})

# Check if biomass data exists, otherwise generate it
if (
  !file.exists(file.path(
    "data",
    "processed",
    "interpolation",
    "weekly_biomasses.csv"
  ))
) {
  system(paste("Rscript", file.path("code", "InterpolateWeekly.R")))
}
# Check if predator selectivity data exists, otherwise generate it.
if (!file.exists(file.path("data", "processed", "predator_selectivity.csv"))) {
  system(paste("Rscript", file.path("code", "CombineMetabarcoding.R")))
}

message("Running ModelForageResponse.R")
# Merge selectivity and biomass ----
# Load and process biomass data

# For the fish prey, use the count data
count_spras <-
  read_csv(
    file.path("data", "raw", "count_spras2022.csv"),
    show_col_types = FALSE
  ) |>
  # Use the bodymass data from BY31 LANDSORTSDJ to estimate the biomass of prey availability based on counts
  left_join(
    read_csv(
      file.path("data", "processed", "interpolation", "weekly_bodymass.csv"),
      show_col_types = FALSE
    ) |>
      filter(station_name == "BY31 LANDSORTSDJ"),
    by = c("node_name", "sample_week")
  ) |>
  na.omit() |>
  mutate(biomass = abundance * bodymass) |>
  select(node_prey = node_name, trawl_id, sample_week, station_name, biomass)

fish_prey <- count_spras |> pull(node_prey) |> unique()
"%!in%" <- Negate("%in%")

# For the zooplankton prey, use the monitoring data
biomass <-
  read_csv(
    file.path("data", "processed", "interpolation", "weekly_biomasses.csv"),
    show_col_types = FALSE
  ) |>
  group_by(sample_week, station_name) |>
  mutate(rel_biomass = biomass / sum(biomass)) |>
  select(node_prey = node_name, sample_week, station_name, biomass) |>
  mutate(trawl_id = NA) |>
  # For the fish prey, use the count data
  filter(node_prey %!in% fish_prey) |>
  rbind(count_spras)

# Load and compute forage ratios
ForageRatios <-
  read_csv(
    file.path("data", "processed", "predator_selectivity.csv"),
    show_col_types = FALSE
  ) |>
  ungroup() |>
  mutate(
    rra_gut = replace_na(rra_gut, 0),
    rra_env = replace_na(rra_env, 0),
    ForageRatio = rra_gut / rra_env
  ) |>
  filter(rra_env > 0) |>
  left_join(
    biomass,
    by = c("sample_week", "station_name", "node_prey", "trawl_id")
  ) |>
  group_by(sample_id) |>
  mutate(rel_biomass = biomass / sum(biomass, na.rm = T)) |>
  ungroup()

# Save the trawl_id for plotting the map later...
ForageRatios |>
  filter(
    !is.na(trawl_id),
    !is.na(biomass),
    node_predator %in% c("Clupea", "Gasterosteus", "Sprattus")
  ) |>
  mutate(
    organism = case_when(
      node_predator == "Gasterosteus" ~ "Gasterosteus aculeatus",
      node_predator == "Clupea" ~ "Clupea harengus",
      node_predator == "Sprattus" ~ "Sprattus sprattus"
    )
  ) |>
  select(organism, trawl_id) |>
  unique() |>
  write_csv(file.path("data", "processed", "trawl_summary.csv"))

# Calculate the average forage ratio for each predator-prey pairs
#average_forage_ratios <-
#  ForageRatios |>
#  filter(!is.na(biomass)) |>
#  group_by(node_predator, node_prey) |>
#  summarise(ForageRatio = mean(ForageRatio, na.rm = T), .groups = "drop")
node_data <- read_csv(
  file.path("data", "raw", "node_data.csv"),
  show_col_types = FALSE
)
fish <- node_data |> filter(type == "fish") |> select(node_name)
zooplankton <- node_data |> filter(type == "zooplankton") |> select(node_name)
phytoplankton <- node_data |>
  filter(type == "phytoplankton") |>
  select(node_name)
neutral_forage_ratios <-
  fish |>
  rename(node_predator = node_name) |>
  cross_join(zooplankton |> rename(node_prey = node_name)) |>
  bind_rows(
    zooplankton |>
      rename(node_predator = node_name) |>
      cross_join(phytoplankton |> rename(node_prey = node_name))
  ) |>
  mutate(ForageRatio = 1)
# Helper functions ----

# Estimate the forage ratio as a density dependant response
# 1: Define a bootstrap function
bootstrap_fit <- function(df, model, coef, prey, predator) {
  df_boot <- df[sample(nrow(df), replace = TRUE), ] # Resample with replacement

  boot_fit <- tryCatch(
    nlsLM(
      model,
      start = list(c = coef[1, 1]),
      data = df_boot,
      lower = c(c_min),
      upper = c(c_max),
      control = nls.lm.control(maxiter = 500)
    ),
    error = function(e) NULL
  )

  if (!is.null(boot_fit)) {
    boot_coef <- summary(boot_fit)$coefficients
    return(tibble(
      node_predator = predator,
      node_prey = prey,
      c = boot_coef[1, 1]
    ))
  } else {
    return(NULL)
  }
}

# 2: Ensure 1000 Valid Bootstrap fits
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
      warning(glue::glue(
        "Reached max attempts for {predator}-{prey}: {length(successful_results)} fits collected."
      ))
      break
    }
  }

  bind_rows(successful_results) |>
    mutate(Iteration = row_number())
}

# 3:Define function to fit model for a single group
fit_group <- function(df) {
  if (nrow(df) <= 3) return(NULL)

  prey <- unique(df$node_prey)
  predator <- unique(df$node_predator)

  # Try all c_values and keep the best fit
  fits <- map_dfr(c_values, function(c_start) {
    fit <- tryCatch(
      nlsLM(
        ForageRatio ~ (1 + c) / (1 + c * rel_biomass),
        start = list(c = c_start),
        data = df,
        lower = c(c_min),
        upper = c(c_max),
        control = nls.lm.control(maxiter = 500)
      ),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      tibble(
        start = c_start,
        residual_dispersion = sum(residuals(fit)^2),
        fit = list(fit)
      )
    } else {
      NULL
    }
  })

  if (nrow(fits) == 0) return(NULL)

  # Select the best fit based on the starting c
  best_fit <- fits |>
    slice_min(residual_dispersion, n = 1) |>
    pull(fit)

  best_fit <- best_fit[[1]]
  coef <- summary(best_fit)$coefficients

  # Plot the residual fit and save as a png
  residuals_df <- tibble(
    fitted = fitted(best_fit),
    residuals = resid(best_fit)
  )

  plot_title <- paste(predator, prey, sep = " - ")
  file_name <- file.path(
    "output",
    "ModelForageResponse",
    "residuals",
    paste0(predator, "_", prey, ".png")
  )

  p <- ggplot(residuals_df, aes(x = fitted, y = residuals)) +
    geom_hline(yintercept = 0, color = "red") +
    geom_point(shape = 21) +
    labs(
      title = plot_title,
      x = "Fitted values",
      y = "Residuals"
    ) +
    theme_bw()

  ggsave(filename = file_name, plot = p, width = 4, height = 4, dpi = 300)

  # Return a tibble, with predator, prey, c and the bootstrapped c
  tibble(
    node_predator = predator,
    node_prey = prey,
    c = coef[1, 1],
    boot = list(run_bootstraps(df, model, coef, prey, predator, n_boot = 1000))
  )
}

# Prepare the output directory
dir.create(
  file.path("output", "ModelForageResponse", "residuals"),
  recursive = T,
  showWarnings = F
)

# Define the model
model <- ForageRatio ~ (1 + c) / (1 + c * rel_biomass)

# Define parameter c ranges and starting values
c_values = seq(-.9, 10, .1)
c_min = -0.99
c_max = 50

# Fit the Model for each predator-prey pair
# Split data by Barcode, and Species_English
df_list <-
  ForageRatios |>
  group_split(node_predator, node_prey, .keep = TRUE)

# Run in parallel with progress bar
plan(multisession)

model_results <-
  future_map(
    df_list,
    fit_group,
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  ) |>
  compact() |>
  bind_rows()

# Save the results ----
model_results |>
  select(-boot) |>
  right_join(neutral_forage_ratios, by = c("node_predator", "node_prey")) |>
  mutate(neutral_forage_ratio = ForageRatio) |>
  select(node_predator, node_prey, neutral_forage_ratio, c) |>
  arrange(node_predator, node_prey) |>
  write_csv(file = file.path("data", "processed", "forage_ratio.csv"))
model_results |>
  select(boot) |>
  unnest(boot) |>
  right_join(
    neutral_forage_ratios |>
      cross_join(tibble(Iteration = 1:1000)),
    by = c("node_predator", "node_prey", "Iteration")
  ) |>
  filter(!is.na(Iteration)) |>
  mutate(neutral_forage_ratio = ForageRatio) |>
  select(node_predator, node_prey, neutral_forage_ratio, c, Iteration) |>
  arrange(node_predator, node_prey) |>
  write_csv(
    file = file.path("data", "processed", "bootstrap_forage_ratio.csv")
  )
# Summarise the confidence intervals from bootstraps
bootstrapped_values <-
  model_results |>
  select(boot) |>
  unnest(boot) |>
  group_by(node_prey, node_predator) |>
  summarise(
    AVG_c = mean(c),
    LOW_c = quantile(c, 0.025),
    HIGH_c = quantile(c, 0.975),
    .groups = "drop"
  ) |>
  mutate(
    c = paste0(
      round(AVG_c, 1),
      "\n[",
      round(LOW_c, 1),
      ";",
      round(HIGH_c, 1),
      "]"
    ),
  ) |>
  select(node_predator, node_prey, c) |> #, h) |>
  arrange(node_predator, node_prey)
if (!dir.exists(file.path("output", "table")))
  dir.create(file.path("output", "table"), recursive = TRUE)
write_csv(
  bootstrapped_values,
  file.path("output", "table", "forage_ratio_parameters.csv")
)
# Quick visualisation -------
rel_biomass_seq <- seq(0, 1, length.out = 100)
boot_prediction <-
  model_results |>
  select(boot) |>
  unnest(boot) |>
  cross_join(tibble(rel_biomass = rel_biomass_seq)) |>
  mutate(
    ForageRatio = (1 + c) / (1 + c * rel_biomass),
    Bgut = ForageRatio * rel_biomass
  ) |>
  group_by(node_predator, node_prey, rel_biomass) |>
  summarise(
    fr_lower = quantile(ForageRatio, 0.025, na.rm = TRUE),
    fr_upper = quantile(ForageRatio, 0.975, na.rm = TRUE),
    Bgut_lower = quantile(Bgut, 0.025, na.rm = TRUE),
    Bgut_upper = quantile(Bgut, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

plot_and_save_curves <- function(p) {
  model_results |>
    select(-boot) |>
    right_join(neutral_forage_ratios, by = c("node_predator", "node_prey")) |>
    cross_join(tibble(rel_biomass = rel_biomass_seq)) |>
    mutate(
      ForageRatio = ifelse(
        !is.na(c),
        (1 + c) / (1 + c * rel_biomass),
        ForageRatio
      )
    ) |>
    filter(
      node_predator == p
    ) |>
    ggplot() +
    geom_hline(yintercept = 1) +
    geom_ribbon(
      data = boot_prediction |>
        filter(
          node_predator == p
        ),
      mapping = aes(ymin = fr_lower, ymax = fr_upper, x = rel_biomass),
      alpha = .2,
      col = "black",
      linetype = 2
    ) +

    geom_point(
      data = ForageRatios |>
        filter(
          node_predator == p
        ) |>
        filter(!is.na(biomass)),
      mapping = aes(x = rel_biomass, y = ForageRatio),
      color = "black",
      shape = 21,
      alpha = .5,
      size = 1
    ) +
    geom_line(
      linewidth = .5,
      col = "red",
      mapping = aes(x = rel_biomass, y = ForageRatio)
    ) +

    facet_wrap(. ~ node_prey, scales = "free", ncol = 5) +

    theme_bw() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 7),
      panel.grid = element_blank(),
      axis.text = element_text(color = "black", size = 7),
      axis.title = element_text(color = "black", size = 8)
    ) +

    scale_x_continuous(breaks = c(0, .5, 1)) +
    labs(x = "Relative Biomass", y = "Forage ratio", title = p)
  name <- paste0("fitted_", p, ".pdf")
  ggsave(
    filename = file.path("output", "ModelForageResponse", name),
    width = 7,
    height = 4
  )
}

predator <- unique(neutral_forage_ratios$node_predator)

walk(predator, ~ plot_and_save_curves(.x))

# Simulation ----
forage_ratio <- function(biomass, c) {
  (1 + c) / (1 + c * biomass)
}
df <- tibble(
  #c = seq(-0.9,10, .01)
  c = c(-.9, -0.7, -0.4, 1, 3, 10)
) |>
  cross_join(tibble(rBe = seq(0, 1, length.out = 200))) |>
  mutate(
    S = forage_ratio(rBe, c),
    rBg = S * rBe
  )

p1 <- df |>
  ggplot(aes(x = rBe, y = S, col = c + 1, group = factor(c))) +
  geom_hline(yintercept = 1) +
  geom_line(linewidth = 1.5) +
  theme_bw() +
  scale_y_log10(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_gradientn(
    colors = c("#fc8d59", "#ffffbf", "#91bfdb"),
    trans = "log10"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_blank(),
    strip.background = element_blank()
  )
p2 <- df |>
  ggplot(aes(x = rBe, y = rBg, col = c + 1, group = factor(c))) +
  geom_abline() +
  geom_line(linewidth = 1.5) +
  theme_bw() +
  scale_color_gradientn(
    colors = c("#fc8d59", "#ffffbf", "#91bfdb"),
    trans = "log10"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_blank(),
    strip.background = element_blank()
  )
plot <- p2 + p1 + plot_layout(guides = "collect", axes = "collect")
ggsave(
  plot = plot,
  filename = file.path("output", "ModelForageResponse", "model.pdf"),
  width = 8,
  height = 4
)
