#!/usr/bin/env Rscript

# List of required packages
required_packages <- c(
  "tidyverse",
  "vegan",
  "ape",
  "rlang",
  "patchwork",
  "tidygraph",
  "ggraph",
  "gridExtra",
  "grid",
  "data.table",
  "spaa",
  "broom",
  "ggrepel",
  "zoo",
  "minpack.lm",
  "furrr",
  "progressr",
  "sf",
  "igraph",
  "fluxweb",
  "abind",
  "renv"
)


# Function to install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing missing package:", pkg))
    install.packages(pkg)
  }
}

# Apply the function to each package
invisible(lapply(required_packages, install_if_missing))
