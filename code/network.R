# Some network thingy
library(igraph)
library(fluxweb)
suppressPackageStartupMessages(library(tidyverse))
library(furrr)
# Function from https://doi.org/10.1111/1365-2656.13447
lw <- function(fluxes, loop = FALSE, parameter = "connectance") {
  #res <- c()
  # The flux matrix
  W.net <- as.matrix(fluxes) #fluxmatrix from fluxweb

  ### Taxon-specific Shannon indices of inflows
  # sum of k species inflows --> colsums
  sum.in <- apply(W.net, 2, sum)

  # Diversity of k species inflows
  # columns divided by the total col sum
  H.in.mat <- t(t(W.net) / sum.in) * t(log(t(W.net) / sum.in))
  H.in.mat[!is.finite(H.in.mat)] <- 0 #converts NaN to 0's
  H.in <- apply(H.in.mat, 2, sum) * -1

  # Effective number of prey or resources = N(R,k)
  # The reciprocal of H(R,k) --> N (R,k) is the equivalent number of prey for species k
  N.res <- ifelse(sum.in == 0, H.in, exp(H.in))

  ### Taxon-specific Shannon indices of outflows
  # sum of k speies outflows --> rowsums
  sum.out <- apply(W.net, 1, sum)

  # Diversity of k species outflows
  # rows divided by the total row sum
  H.out.mat <- (W.net / sum.out) * log(W.net / sum.out)
  H.out.mat[!is.finite(H.out.mat)] <- 0 #converts NaN to 0's
  H.out <- apply(H.out.mat, 1, sum) * -1

  # Effective number of predators or consumers = N(C,k)
  # The reciprocal of H(C,k) --> N (C,k) is the equivalent number of predators for species k
  N.con <- ifelse(sum.out == 0, H.out, exp(H.out))

  ### Quantitative Weighted connectance
  no.species <- ncol(W.net)

  # The weighted link density (LDw) is:
  # In the weighted version the effective number of predators for species i is weighted by i's
  # contribution to the total outflow the same is the case for the inflows
  tot.mat <- sum(W.net)
  # LD.w <- (sum((sum.in/tot.mat)*N.res) + sum((sum.out/tot.mat)*N.con))/2
  # equivalent to next formula, but next one is closer to manuscript
  LD <- 1 / (2 * tot.mat) * (sum(sum.in * N.res) + sum(sum.out * N.con))

  #Weighted connectance
  lwC <- LD / ifelse(loop, no.species, no.species - 1)

  # positional.index
  #pos.ind<- sum.in*N.res/(sum.in*N.res+sum.out*N.con) #postional index
  #basal.sp<- pos.ind[pos.ind==0] #basal species = 0
  #top.sp<- pos.ind[pos.ind==1] #defintion according to Bersier et al. 2002 top species = [0.99, 1]

  #con.sp<-length(pos.ind)-length(basal.sp)# all consumer taxa except basal
  # weighted quantitative Generality
  lwG <- sum(sum.in * N.res / sum(W.net))

  #res.sp<- length(pos.ind)-length(top.sp)
  # weighted quantitative Vulnerability
  lwV <- sum(sum.out * N.con / sum(W.net))

  if (parameter == "connectance") return(lwC)
  if (parameter == "generality") return(lwG)
  if (parameter == "vulnerability") return(lwV)
}
source("./code/HelpMyAnalyses.R")
timeseries_fluxes <- read_csv(
  file = file.path("data", "analyses", "timeseries_fluxes.csv"),
  show_col_types = FALSE
)
timeseries_null <- read_csv(
  file = file.path("data", "analyses", "timeseries_null.csv"),
  show_col_types = FALSE
)
timeseries_null |> str()
timeseries_fluxes |>
  filter(sample_week == "2008-01-07") |>
  select(predator, prey, mean) |>
  pivot_wider(values_from = mean, names_from = prey, values_fill = 0) |>
  column_to_rownames("predator") |>
  as.matrix() |>
  t()

connectance_df <-
  timeseries_fluxes |>
  filter(isoweek(sample_week) %in% 2:51, year(sample_week) %in% 2008:2023) |>

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

connectance_df |>
  pivot_longer(2:4, values_to = "values", names_to = "parameter") |>
  group_by(iso_week, parameter) |>
  summarise(avg = mean(values), SD = sd(values), .groups = "drop") |>
  ggplot(aes(x = iso_week, y = avg, ymin = avg - SD, ymax = avg + SD)) +
  geom_ribbon(alpha = .2) +
  geom_line() +
  facet_grid(parameter ~ ., scales = "free")

connectance_df |>
  filter(season != "Winter") |>
  pivot_longer(2:4, values_to = "values", names_to = "parameter") |>
  group_by(year = year(sample_week), parameter) |>

  summarise(value = mean(values), .groups = "drop") |>
  ggplot(aes(x = year, y = value)) +
  geom_point() +
  facet_grid(parameter ~ ., scales = "free") +
  geom_smooth()

connectance_null <-
  timeseries_null |>
  filter(isoweek(sample_week) %in% 2:51, year(sample_week) %in% 2008:2023) |>

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

connectance_df |>
  pivot_longer(2:4, values_to = "values", names_to = "parameter") |>
  mutate(selectivity = T) |>
  rbind(
    connectance_null |>
      pivot_longer(2:4, values_to = "values", names_to = "parameter") |>
      mutate(selectivity = F)
  ) |>
  group_by(iso_week, parameter, selectivity) |>
  summarise(avg = mean(values), SD = sd(values), .groups = "drop") |>
  mutate(parameter = paste("lw", parameter)) |> 
  ggplot(aes(
    x = iso_week,
    y = avg,
    ymin = avg - SD,
    ymax = avg + SD,
    col = selectivity,
    fill = selectivity
  )) +
  geom_ribbon(alpha = .2, col = NA) +
  geom_line() +
  facet_grid(parameter ~ ., scales = "free", switch = "y") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text = element_text(color = "black")) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(
    breaks = seq(1, 52.1775, 4.348125),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = c("black", "#B80C09")) +
  scale_color_manual(values = c("black", "#B80C09"))
ggsave("./output/figure/network.pdf", width = 5, height = 6)
