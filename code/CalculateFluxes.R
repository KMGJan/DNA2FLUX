#!/usr/bin/env Rscript


# Define Functions ----

#' Filter data by sample week and station
#' 
getStationDate <- function(data, date, station) {
  
  require(tidyverse)
 
   data |> 
    filter(sample_week == date(date),
           station_name == station)
}



#' Temperature metabolic constant
#' For specific date and station
#' 
getTempKonstant <- function(data, date, station) {
  
  require(tidyverse)
  
  temp <- temperature |> 
    getStationDate(date, station) |> 
    pull(temperature)
  boltz <- 0.00008617343 # Boltzmann constant
  tkonst <- 0.69/(boltz*(273.15+temp)) #Temperature metabolic constant
  return(tkonst)
}
#' Example
#'    getTempKonstant(temperature, "2009-02-16", "BY31 LANDSORTSDJ")



#' Get Node Data
#' For specific date and station
#' Calculate Metabolic loss
#' 
getNodeData <- function(node_data, weekly_biomasses,
                        weekly_bodymass, temperature,
                        date, station) {
  
  require(tidyverse)
  
  node_data |> 
    select(node_name, efficiencies, intercept, slope) |> 
    # Add biomasses
    left_join(select(getStationDate(weekly_biomasses, date, station),
                     node_name, biomass),
              by = join_by(node_name)) |>
    mutate(biomass = replace_na(biomass, 0)) |>
    # Add bodymasses
    left_join(select(getStationDate(weekly_bodymass, date, station),
                     node_name, bodymass),
              by = join_by(node_name)) |>
    mutate(bodymass = replace_na(bodymass, 1)) |> 
    # Calculate metabolic losses
    mutate(losses = exp(slope * log(bodymass) + intercept - getTempKonstant(temperature, date, station)),
           losses = ifelse(is.infinite(losses), 0, losses))
  
}
#' Example
#' getNodeData(node_data, weekly_biomasses,
#'            weekly_bodymass, temperature,
#'            "2009-02-16", "BY31 LANDSORTSDJ") 






tidyFluxing <- function(graph) {
  
  require(fluxweb)
  require(igraph)
  require(tidyverse)
  
  fluxes <- 
    fluxing(mat = as_adjacency_matrix(graph, attr = "weight", sparse = FALSE),
            biomasses = pull(graph, biomass),
            losses = pull(graph, losses),
            efficiencies = pull(graph, efficiencies),
            bioms.prefs = FALSE,
            ef.level = "pred",
            bioms.losses = TRUE)
}



#' Wrapper
#' Merge all data and calculate fluxes
dna2flux <- function(forage_ratio, node_data,
                     weekly_biomasses, weekly_bodymass,
                     temperature, date, station) {
  
  require(tidyverse)
  require(tidygraph)
  
  forage_ratio |> 
    filter(is.na(node_predator) == F) |> 
    left_join(select(getStationDate(weekly_biomasses, date, station),
                     node_name, biomass),
              by = join_by(node_prey == node_name)) |> 
    group_by(node_predator) |> 
    mutate(rel_biomass = biomass / sum(biomass, na.rm = T)) |> 
    mutate(forage_ratio =ifelse(!is.na(a) & !is.na(h),
                                (a * rel_biomass) / (1 + a * h * rel_biomass) / (rel_biomass),
                                average_forage_ratio),
           forage_ratio = ifelse(is.na(forage_ratio) == T,
                                 0,
                                 forage_ratio),
           weight = (rel_biomass * forage_ratio) / sum(rel_biomass * forage_ratio, na.rm = T)) |> 
    ungroup() |> 
    
    # Make Table Graph
    select(node_predator, node_prey, forage_ratio, weight) |> 
    as_tbl_graph() |>
    activate(edges) |> 
    mutate(weight = replace_na(weight, 0)) |> 
    
    # Add node data  
    activate(nodes) |> 
    left_join(getNodeData(node_data, weekly_biomasses,
                          weekly_bodymass, temperature,
                          date, station),
              by = join_by(name == node_name)) |> 
    tidyFluxing() |> 
    as_tbl_graph() |> 
    activate(nodes) |> 
    left_join(getNodeData(node_data, weekly_biomasses,
                          weekly_bodymass, temperature,
                          date, station),
              by = join_by(name == node_name)) 
}
#' Example
#' dna2flux(forage_ratio, node_data,
#'          weekly_biomasses, weekly_bodymass,
#'          temperature, "2009-02-16", "BY31 LANDSORTSDJ")







# Execute functions --------------------------------------------------

# Load data
library(tidyverse)
library(ggraph)
temperature <- read_csv(file = file.path("data", "Processed", "interpolation", "temperature.csv"))
weekly_biomasses <- read_csv(file = file.path("data", "processed", "interpolation", "weekly_biomasses.csv"))
weekly_bodymass <- read_csv(file = file.path("data", "processed", "interpolation", "weekly_bodymass.csv"))
node_data <- read_csv(file = file.path("data", "raw", "node_data.csv"))
forage_ratio <- read_csv(file = file.path("data", "processed", "forage_ratio.csv"))
bootstrap_forage_ratio <- read_csv(file = file.path("data", "processed", "bootstrap_forage_ratio.csv"))

# Run
dna2flux(forage_ratio, node_data, weekly_biomasses, weekly_bodymass,
         temperature, "2009-02-16", "BY31 LANDSORTSDJ") |> 
  ggraph() +
  geom_edge_link(aes(width = weight), arrow = arrow(length = unit(3, 'mm')),
                 alpha = 0.2) +
  geom_node_point(aes(size = biomass)) +
  geom_node_text(aes(label = name), angle = -90, size = 3, nudge_y = -0.2) +
  theme_graph() +
  #scale_color_manual(values = colors) +
  theme(legend.position = "none")






