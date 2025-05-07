#' Filter data by sample week and station
#' 
#' This function filters the given data for a specific station and sample week. It ensures that the date provided is within the same week (start of the week) for the calculation of fluxes.
#'
#' @param data A data frame containing the data with columns `sample_week` and `station_name`.
#' @param date A Date object representing the specific date within the week to calculate fluxes.
#' @param station A character string representing the station name for which to filter the data.
#' 
#' @return A data frame filtered by the specified date (adjusted to the start of the week) and station name.
#' temperature <- data.frame(temperature = c(20, 10), sample_week = c(as.Date("2022-06-01"), as.Date("2022-06-08")), station_name = rep("BY31 LANDSORTSDJ", 2))
#' getStationDate(data = temperature, date = as.Date("2022-06-01"), station = "BY31 LANDSORTSDJ")
#' 
getStationDate <- function(data, date, station) {
  # Step to ensure that any date can be used within a week to calculate fluxes from the week
  sample_date <- floor_date(date(date), unit = "week", week_start = 1)
  
  data |> 
    filter(sample_week == date(date),
           station_name == station)
}

#' Calculate the temperature metabolic constant
#' 
#' This function calculates the temperature metabolic constant for a specific date and station based on temperature data. The constant is calculated using the Boltzmann constant and the temperature in Celsius.
#' 
#' @param data A data frame containing temperature data with columns `sample_week`, `station_name`, and `temperature`.
#' @param date A Date object representing the specific date to extract temperature data.
#' @param station A character string representing the station name for which to extract temperature data.
#' 
#' @return A numeric value representing the temperature metabolic constant (K_T) for the given date and station.
#' @examples
#' temperature <- data.frame(temperature = c(20, 10), sample_week = c(as.Date("2022-06-01"), as.Date("2022-06-08")), station_name = rep("BY31 LANDSORTSDJ", 2))
#' getTempKonstant(data = temperature, date = as.Date("2022-06-01"), station = "BY31 LANDSORTSDJ")
#' 
getTempKonstant <- function(data, date, station) {
  # Using the function getStationDate, get the temperature at a certain date and station
  temp <- temperature |> 
    getStationDate(date, station) |> 
    pull(temperature)
  
  # Boltzmann constant
  boltz <- 0.00008617343
  
  #Temperature metabolic constant
  tkonst <- 0.69 / (boltz * (273.15 + temp)) 
  return(tkonst)
}

#' Get Node Metadata with Dynamic Biomass, Bodymass, and Metabolic Losses
#'
#' This function returns a complete node data table with dynamic biomass, bodymass, and estimated metabolic losses based on temperature for a given sampling date and station. 
#'
#' @param node_data A data frame with node metadata. Must include columns: `node_name`, `efficiencies`, `intercept`, `slope`, `trophic_level`, `horizontal_position`, and `color`.
#' @param weekly_biomasses A data frame containing biomass values per node per week and station.
#' @param weekly_bodymass A data frame containing bodymass values per node per week and station.
#' @param temperature A data frame containing temperature values per week and station.
#' @param date A `Date` object representing the current sampling date.
#' @param station A character string indicating the station name.
#'
#' @return A data frame that includes node metadata merged with week- and station-specific biomass, bodymass, and computed metabolic losses.
#'
#' @examples
#' getNodeData(node_data, weekly_biomasses, weekly_bodymass, temperature, date = as.Date("2022-06-01"), station = "BY31 LANDSORTSDJ")
#'
getNodeData <- function(node_data, weekly_biomasses, weekly_bodymass, temperature, date, station) {

  node_data |>
    select(node_name, efficiencies, intercept, slope, energy_density_ww, trophic_level, horizontal_position, color) |> 
    
    # Using getStationDate, join the organisms biomass at a given time and location
    left_join(select(getStationDate(data = weekly_biomasses,
                                    date = date,
                                    station = station),
                     node_name, biomass),
              by = join_by(node_name)) |>
    
    # If the biomass is missing, it means that it equals to 0
    mutate(biomass = replace_na(biomass, 0)) |>
    
    # Using getStationDate, join the organisms bodymass at a given time and location
    left_join(select(getStationDate(data = weekly_bodymass,
                                    date = date,
                                    station = station),
                     node_name, bodymass),
              by = join_by(node_name)) |>
    
    # This avoids to have NA, but won't impact caclulations later on
    mutate(bodymass = replace_na(bodymass, 1)) |> 

    # Calculate temperature-corrected metabolic losses
    mutate(losses = exp(slope * log(bodymass) + intercept - getTempKonstant(data = temperature,
                                                                            date = date,
                                                                            station = station)),
           losses = ifelse(is.infinite(losses), 0, losses))
}

#' Compute Trophic Fluxes Using `fluxweb::fluxing` from a `tbl_graph` Object
#'
#' This function calculates energy fluxes in a food web using `fluxweb::fluxing()`. It extracts the adjacency matrix and node attributes from a `tbl_graph` object and ensures that node-level vectors (e.g. biomass, losses) are correctly aligned with the matrix node order.
#'
#' @param graph A `tbl_graph` object (from the `tidygraph` package), where:
#'   - Edge weights represent interaction strengths (Wij) from prey to predators.
#'   - Node attributes must include `biomass`, `losses`, and `efficiencies`, with names matching the nodes in the graph.
#'
#' @return A numeric matrix of estimated trophic fluxes (predators in rows, prey in columns).
#'
#' @details
#' This function acts as a tidy wrapper around `fluxweb::fluxing()`, and:
#' - Extracts the weighted adjacency matrix with `tidygraph::as_adjacency_matrix()`, then transposes it so predators are in rows.
#' - Extracts node attributes and orders them to match the matrix's node order (`igraph::V(graph)$name`).
#' - Calls `fluxweb::fluxing()` with the following options:
#'   - `bioms.prefs = FALSE`: no biomass preference in foraging.
#'   - `ef.level = "pred"`: efficiencies apply to predators.
#'   - `bioms.losses = TRUE`: metabolic losses scale with biomass.
#' - The resulting flux matrix is transposed again so that the output has predators in rows and prey in columns.
#'
#' @importFrom fluxweb fluxing
#' @importFrom tidygraph as_adjacency_matrix
#' @importFrom igraph V
#' @importFrom dplyr pull filter arrange match
#'
#' @seealso [fluxweb::fluxing()], [tidygraph::as_adjacency_matrix()], [igraph::V()]
#'
#' @examples
#' \dontrun{
#'   tidyFluxing(graph)
#' }
#' 
tidyFluxing <- function(graph, population_growth = TRUE) {
  node_names <- V(graph)$name
  node_data <- as_tibble(graph, what = "vertices")
  
  # Ensure correct order of attributes
  node_data_ordered <- node_data |> 
    filter(name %in% node_names) |> 
    arrange(match(name, node_names))
  
  # Build adjacency matrix (prey x predator)
  adj_mat <- t(as_adjacency_matrix(graph, attr = "weight", sparse = FALSE))
  
  # Choose the loss vector depending on population growth
  losses <- if (population_growth) {
    node_data_ordered$population_losses
  } else {
    node_data_ordered$losses
  }
  
  # Compute fluxes and transpose to predator x prey
  flux_mat <- fluxing(
    mat = adj_mat,
    biomasses = node_data_ordered$biomass,
    losses = losses,
    efficiencies = node_data_ordered$efficiencies,
    bioms.prefs = FALSE,
    ef.level = "pred",
    bioms.losses = !population_growth
  )
  
  return(t(flux_mat))
}

#' Calculate Trophic Fluxes from Forage Ratios
#'
#' This function estimates energy or mass fluxes between nodes in a food web based on predator-prey forage ratios, prey biomasses, and metabolic losses. It returns either a numeric flux matrix or a `tbl_graph` with flux values as edge weights.
#'
#' @param forage_ratio A data frame containing forage ratio parameters, with columns:
#'   - `node_predator`, `node_prey`: predator-prey node identifiers,
#'   - `a`, `h`, `average_forage_ratio`: parameters for functional response.
#' @param node_data A node-level metadata table including model parameters such as trophic level, etc.
#' @param weekly_biomasses A data frame of biomasses per node, sample week, and station.
#' @param weekly_bodymass A data frame of bodymass values per node, sample week, and station.
#' @param temperature A data frame with temperature values per week and station.
#' @param date The sample week (date) for which fluxes should be estimated.
#' @param station The name of the station/site for filtering biomass and temperature data.
#' @param as_graph Logical. If `TRUE`, returns a `tbl_graph` with fluxes on edges; otherwise returns a matrix of fluxes in `kJ/day/m²`.
#'
#' @return A numeric flux matrix (if `as_graph = FALSE`) or a `tidygraph::tbl_graph` with nodes and weighted edges (if `as_graph = TRUE`).
#'
#' @details
#' The function:
#' 1. Computes relative prey biomasses per predator.
#' 2. Calculates forage ratios using a Type II functional response.
#' 3. Converts forage ratios into normalized weights for each predator.
#' 4. Builds a `tbl_graph` with weighted interactions and node attributes.
#' 5. Applies `tidyFluxing()` to compute fluxes based on biomass, losses, and efficiencies.
#' 6. Converts fluxes from J/s/m² to kJ/day/m² using a constant (86.4).
#'
#'
#' @examples
#' \dontrun{
#'   flux_matrix <- dna2flux(forage_ratio, node_data, biomasses, bodymass,
#'                           temperature, date = "2010-06-01", station = "BY31 LANDSORTSDJ")
#'   flux_graph <- dna2flux(forage_ratio, node_data, biomasses, bodymass,
#'                          temperature, date = "2010-06-08", station = "BY31 LANDSORTSDJ",
#'                          as_graph = TRUE)
#' }
#' 
dna2flux <- function(forage_ratio, node_data, weekly_biomasses, weekly_bodymass, temperature, date, station,  as_graph = FALSE, presence_absence = FALSE, population_growth = TRUE) {
  
    # Change this to: Xi = Pi*Ri + ((Delta_B/Delta_t)*ED)
    if(population_growth == TRUE){
      # Change this to: Xi = Pi*Ri + ((Delta_B/Delta_t)*ED)
      #Delta_B = Biomass week n+1 - Biomass week
      date_next_week = date + 7
      
      node_values_next_week <- getNodeData(node_data = node_data,
                                           weekly_biomasses = weekly_biomasses,
                                           weekly_bodymass = weekly_bodymass,
                                           temperature = temperature,
                                           date = date_next_week,
                                           station = station) |>
        select(node_name, "biomass_next_week" = biomass)
      
      node_values_this_week <- getNodeData(node_data = node_data,
                                           weekly_biomasses = weekly_biomasses,
                                           weekly_bodymass = weekly_bodymass,
                                           temperature = temperature,
                                           date = date,
                                           station = station)
        
      node_values <-  node_values_this_week |>
        left_join(node_values_next_week, by = join_by(node_name)) |> 
        mutate(
          individual_rate = losses,
          dP = biomass_next_week - biomass,  
          dt = 604800, # 7 days in seconds
          basal_metabolism = biomass * individual_rate,
          recruitment = pmax(0, dP/dt * energy_density_ww), # This results in positive values or zero
          mortality = pmin(0, dP * individual_rate), # This results in negative values or zero
          population_losses = basal_metabolism + recruitment + mortality
        )
        
    } else {
    node_values <- getNodeData(node_data = node_data,
                               weekly_biomasses = weekly_biomasses,
                               weekly_bodymass = weekly_bodymass,
                               temperature = temperature,
                               date = date,
                               station = station)
  }
  
  tbl <- 
     forage_ratio |> 
     filter(!is.na(node_predator)) |> 
     left_join(
       select(getStationDate(data = weekly_biomasses,
                             date = date,
                             station = station),
              node_name, biomass),
       by = join_by(node_prey == node_name)
     ) |> 
     group_by(node_predator) |> 
     mutate(rel_biomass = biomass / sum(biomass, na.rm = TRUE),
     # Calculate weight depending on presence_absence
       weight = if (presence_absence) {
         presence <- if_else(average_forage_ratio > 0, 1, 0)
         (rel_biomass * presence) / sum(rel_biomass * presence, na.rm = TRUE)
       } else {
         forage_ratio <- case_when(
           !is.na(a) & !is.na(h) ~ (a * rel_biomass) / (1 + a * h * rel_biomass) / rel_biomass,
           TRUE ~ average_forage_ratio
         )
         forage_ratio <- replace_na(forage_ratio, 0)
         (rel_biomass * forage_ratio) / sum(rel_biomass * forage_ratio, na.rm = TRUE)
       }
     ) |> ungroup() |> 
   
     # Make Table Graph
     select(node_predator, node_prey, weight) |> 
     as_tbl_graph() |>
     activate(edges) |> 
     mutate(weight = replace_na(weight, 0)) |> 
     
     # Add node data  
     activate(nodes) |> 
     left_join(
       node_values,
       by = join_by(name == node_name))
   
   mat <- tbl |> 
     tidyFluxing(population_growth = population_growth) * 86.4 # From J/second/m2 to kJ/day/m2 

  graph <- 
    mat |>
    as_tbl_graph() |> 
    activate(nodes) |>
    left_join(
      getNodeData(node_data = node_data,
                  weekly_biomasses = weekly_biomasses,
                  weekly_bodymass = weekly_bodymass,
                  temperature = temperature,
                  date = date,
                  station = station),
              by = join_by(name == node_name))
    
  if (as_graph == TRUE) {
    return(graph)
  } else {
    return(mat)
  }
}


#' Bootstrap Trophic Fluxes
#'
#' Applies `dna2flux()` to each bootstrap replicate of the forage ratio input to estimate the variability in trophic fluxes. Returns a 3D array with flux estimates per bootstrap iteration.
#'
#' @param bootstrap_forage_ratio A data frame of forage ratio bootstrap replicates. Must contain an `Iteration` column to distinguish bootstrap samples.
#' @param node_data Metadata about each node, including model parameters such as slope, intercept, trophic level, efficiencies, etc.
#' @param weekly_biomasses A data frame with weekly biomass per node and station.
#' @param weekly_bodymass A data frame with bodymass per node and station.
#' @param temperature A data frame with weekly temperature per station.
#' @param date The sample date (or week) for which fluxes are calculated.
#' @param station The station name to filter relevant environmental and biological data.
#' @param as_graph Logical. Passed to `dna2flux()`. If `TRUE`, returns a list of `tbl_graph` objects. Typically `FALSE` for bootstrapping to return numeric matrices.
#'
#' @return A 3D array (`n_predators` x `n_prey` x `n_iterations`) containing trophic fluxes in kJ/day/m². Each slice along the third dimension corresponds to one bootstrap iteration.
#'
#' @details
#' - Uses `purrr::possibly()` to safely call `dna2flux()` in case of failure for individual replicates.
#' - Drops any NULL outputs returned due to failure in `dna2flux()` via `purrr::keep()`.
#' - Combines valid matrices into a single 3D array using `abind::abind()`.
#'
#' @examples
#' \dontrun{
#' boot_array <- bootstrapFluxes(bootstrap_forage_ratio, node_data, biomasses, bodymass, temperature, date = "2012-06-01", station = "BY31 LANDSORTSDJ")
#' }
#'
bootstrapFluxes <- function(bootstrap_forage_ratio, node_data, weekly_biomasses, weekly_bodymass, temperature, date, station, as_graph = FALSE, presence_absence = FALSE, population_growth = TRUE) {
  
  safe_dna2flux <- possibly(dna2flux, otherwise = matrix(nrow = 24, ncol = 24))
  
  bootstrap_forage_ratio |>
    group_by(Iteration) |> 
    group_split() |> 
    map(function(group) {
      safe_dna2flux(forage_ratio = group,
                    node_data = node_data,
                    weekly_biomasses = weekly_biomasses,
                    weekly_bodymass = weekly_bodymass,
                    temperature = temperature,
                    date = date,
                    station = station,
                    as_graph = FALSE,
                    presence_absence = FALSE,
                    population_growth = population_growth)
    }) |> 
    keep(~ !is.null(.)) |> 
    abind::abind(along = 3)
}

#' Cache bootstrapped trophic flux arrays
#'
#' This function checks if a bootstrapped trophic flux array has already been cached for a given station and date. If not, it computes the flux array and saves it as an `.rds` file to avoid re-computation in future runs.
#'
#' @param cache.dir Character. Path to the directory where cached `.rds` files are stored.
#' @param bootstrap_forage_ratio A bootstrapped forage ratio object used for flux estimation.
#' @param node_data A data frame containing node-specific attributes (e.g., metabolic rates).
#' @param weekly_biomasses A data frame or list with weekly biomass values for each node.
#' @param weekly_bodymass A data frame or list with weekly body mass values for each node.
#' @param temperature Numeric. The temperature (°C) associated with the trophic flux computation.
#' @param date Character. The date corresponding to the flux estimation (e.g., `"2022-10-15"`).
#' @param station Character. The station identifier corresponding to the sampling location.
#' @param as_graph Logical. Whether the resulting flux array should be returned as a `tbl_graph` object (default `FALSE`). **Note**: this argument is passed to `bootstrapFluxes()` but is not used directly here.
#'
#' @return No return value. This function is called for its side effect of writing a `.rds` file to `cache.dir` if the file does not already exist.
#'
cacheMyFluxes <- function(cache.dir, bootstrap_forage_ratio, node_data, weekly_biomasses, weekly_bodymass, temperature, date, station, as_graph = FALSE, presence_absence = FALSE, population_growth = TRUE) {
  
  cache_file <- file.path(cache.dir, paste0("flux_", station, "_", date, ".rds"))
  if (!dir.exists(cache.dir)) dir.create(cache.dir)

  if (!file.exists(cache_file)) {
    bootstrapFluxes(bootstrap_forage_ratio = bootstrap_forage_ratio,
                    node_data = node_data,
                    weekly_biomasses = weekly_biomasses,
                    weekly_bodymass = weekly_bodymass,
                    temperature = temperature,
                    date = date,
                    station = station,
                    as_graph = as_graph,
                    presence_absence = presence_absence,
                    population_growth = population_growth) |> 
      write_rds(cache_file)
  }
}

#' Estimate Weekly Trophic Fluxes with Bootstrapped Confidence Intervals
#'
#' This function computes weekly trophic fluxes by applying statistical summaries (mean, 2.5%, and 97.5% quantiles) across bootstrapped flux estimates. Optionally uses caching to avoid redundant calculations.
#'
#' @param bootstrap_forage_ratio A data frame of bootstrap replicates for forage ratios. Must include an `Iteration` column to distinguish samples.
#' @param node_data Metadata for each node, including parameters such as metabolic slope, intercept, trophic level, efficiencies, etc.
#' @param weekly_biomasses Weekly biomass data by node and station.
#' @param weekly_bodymass Weekly body mass data by node and station.
#' @param temperature Weekly temperature data by station.
#' @param date The sampling date (or start of week) for which fluxes are estimated.
#' @param station The station name used to filter environmental and biological data.
#' @param cache.dir Optional path to a directory for caching bootstrapped flux matrices. If `FALSE` (default), no caching is used.
#'
#' @return A `tbl_graph` object representing the flux network with nodes enriched with ecological and environmental metadata. Edges contain mean fluxes and 95% confidence intervals.
#'
#' @details
#' - Internally calls `bootstrapFluxes()` and summarizes the resulting 3D array.
#' - Confidence bounds are computed as 2.5% and 97.5% quantiles of bootstrapped fluxes.
#' - If `cache.dir` is specified, results are stored and retrieved as RDS files named by station and date.
#' - Resulting graph is ready for further analysis or visualization using the `tidygraph` and `ggraph` ecosystems.
#'
#'
fluxingWithConfidence <- function(bootstrap_forage_ratio, node_data, weekly_biomasses, weekly_bodymass, temperature, date, station, cache.dir = F, presence_absence = FALSE, population_growth = TRUE, ...) {
  
  # read boostrap_list from cache if it exists
  if (cache.dir != FALSE) {
    cacheMyFluxes(cache.dir, bootstrap_forage_ratio, node_data, weekly_biomasses, weekly_bodymass, temperature, date, station, as_graph = FALSE, presence_absence = FALSE, population_growth = TRUE)
    cache_file <- file.path(cache.dir, paste0("flux_", station, "_", date, ".rds"))
    bootstrap_array <- read_rds(cache_file)
  } else {
    bootstrap_array <- bootstrapFluxes(bootstrap_forage_ratio = bootstrap_forage_ratio,
                                      node_data = node_data,
                                      weekly_biomasses = weekly_biomasses,
                                      weekly_bodymass = weekly_bodymass,
                                      temperature = temperature,
                                      date = date,
                                      station = station,
                                      as_graph = as_graph,
                                      presence_absence = presence_absence,
                                      population_growth = population_growth)
  }
  
  
  summaries <- list(
    flux_mean = as_function(~mean(.x, na.rm = TRUE)),
    flux_lower = as_function(~quantile(.x, probs = 0.025, na.rm = TRUE)),
    flux_upper = as_function(~quantile(.x, probs = 0.975, na.rm = TRUE))
  )
  

  flux_long <- imap(summaries, function(f, name) {
    apply(bootstrap_array, c(1, 2), f) |>
      as.data.frame() |>
      rownames_to_column("predator") |>
      pivot_longer(-predator, names_to = "prey", values_to = name)
  }) |>
    reduce(left_join, by = c("predator", "prey"))
  
  flux_long |>
    as_tbl_graph() |> 
    activate(nodes) |> 
    left_join(
      getNodeData(node_data, weekly_biomasses,
                  weekly_bodymass, temperature,
                  date, station),
      by = join_by(name == node_name)
    )
}

#' Extract and reshape trophic flux matrices from a graph
#'
#' This function extracts flux matrices (mean, upper, lower) stored as edge attributes in a `tbl_graph` object, reshapes them into long format, and combines them into a single data frame. It also adds the corresponding sampling week to each flux record.
#'
#' @param graph A `tbl_graph` object where nodes have a `sample_week` attribute  and edges have `flux_mean`, `flux_upper`, and `flux_lower` attributes.
#'
#' @return A tibble in long format with columns: `predator`, `prey`, `flux_mean`, `flux_upper`,  `flux_lower`, and `sample_week`.
#'
#' @details The function internally converts adjacency matrices (for mean, upper, and lower fluxes) to long format using a helper function. Only non-NA flux values are kept.
#'
extract_flux_long <- function(graph) {
  week <- graph |>
    activate(nodes) |>
    as_tibble() |>
    pull(sample_week) |>
    unique()

  # Helper function to convert matrix to long format
  adj_to_long <- function(graph, attr, value_name) {
    as.matrix(as_adj(graph, attr = attr, sparse = FALSE)) |>
      as.data.frame() |>
      rownames_to_column("predator") |>
      pivot_longer(-predator, names_to = "prey", values_to = value_name)
  }
  # Get all three matrices as long-form data frames
  adj_long_mean <- adj_to_long(graph, "flux_mean", "flux_mean")
  adj_long_upper <- adj_to_long(graph, "flux_upper", "flux_upper")
  adj_long_lower <- adj_to_long(graph, "flux_lower", "flux_lower")
  
  # Join and annotate
  flux_long <- adj_long_mean |>
    filter(!is.na(flux_mean)) |>
    left_join(adj_long_upper, by = c("predator", "prey")) |>
    left_join(adj_long_lower, by = c("predator", "prey")) |>
    mutate(sample_week = week)
  
  return(flux_long)
}

#' Aggregate trophic flux arrays over time
#'
#' This function reads multiple 3D trophic flux arrays (predator x prey x iteration) from `.rds` files, stacks them into a 4D array (predator x prey x iteration x time), averages fluxes across the time dimension, and then summarises the result across bootstrap iterations to compute mean fluxes and 95% confidence intervals.
#'
#' @param cache.dir A character vector of file paths to `.rds` files containing 3D flux arrays (predator x prey x iteration).
#'
#' @return A tibble with columns `predator`, `prey`, `mean`, `lower`, and `upper`, summarizing the average flux and its 95% confidence interval across iterations.
#'
#' @details
#' - The function assumes that each `.rds` file contains a 3D array.
#' - Arrays are combined into a 4D structure to allow aggregation over time (weeks/months).
#' - The function internally uses `abind::abind()` for stacking and `purrr::map()` for reading files.
#'
aggregateFluxes <- function(cache.dir) {
  
  # Read the RDS files into a list of 3D arrays and combine into one 4D array (predator x prey x iteration x time)
  flux_array <- map(cache.dir, read_rds) |> 
    abind(along = 4)

  # Average over time (dimension 4)
  flux_mean_aggregate <- apply(flux_array, c(1, 2, 3), mean, na.rm = TRUE)
  
  # Turn the 3D array into a tidy dataframe
  flux_df <- as.data.frame.table(flux_mean_aggregate, responseName = "flux") |>
    rename(predator = Var1, prey = Var2, iteration = Var3) |>
    mutate(
      flux = as.numeric(flux),
      iteration = as.integer(iteration)
    )
  
  # Summarize across iterations (confidence intervals)
  aggregated_fluxes <-
    flux_df |>
    group_by(predator, prey) |>
    summarise(
      mean = mean(flux, na.rm = TRUE),
      lower = quantile(flux, 0.025, na.rm = TRUE),
      upper = quantile(flux, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(aggregated_fluxes)
}


