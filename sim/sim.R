# Simulations like in Tokar 2018
# Arseniy Khvorov
# Created 2019/12/13
# Last edit 2019/12/13

library(tidyverse)
library(Rcpp)
library(ggdark) # devtools::install_github("khvorov45/ggdark")

# Directories to be used later
sim_dir <- "sim"

sourceCpp(file.path(sim_dir, "sim.cpp"))

# Functions ===================================================================

# Finds density coefficient to result in a particular overall proportion
find_dens_coef <- function(mean, sd, overall_prop, npop, n_days) {
  overall_prop * npop / sum(dnorm(1:n_days, mean, sd))
}

# Counts with a particular overall proportion
get_counts <- function(mean, sd, overall_prop, npop, n_days) {
  dnorm(1:n_days, mean, sd) *
    find_dens_coef(mean, sd, overall_prop, npop, n_days)
}

# Plots a distribution
plot_distr <- function(mean, sd, overall_prop, npop, n_days) {
  tibble(
    day = 1:n_days,
    counts = get_counts(mean, sd, overall_prop, npop, n_days)
  ) %>%
    ggplot(aes(day, counts)) +
    dark_theme_bw(verbose = FALSE) +
    geom_line()
}

# Simulates a population
sim_pop <- function(n_days, init_pop_size, nvac, nflu_novac, ve, lag) {
  sim_pop_cpp(n_days, init_pop_size, nvac, nflu_novac, ve, lag) %>%
    as_tibble()
}

# Script ======================================================================
sourceCpp(file.path(sim_dir, "sim.cpp"))

nsam <- 1e6
n_days <- 300
pop <- sim_pop(
  n_days = n_days,
  init_pop_size = nsam,
  nvac = get_counts(100, 50, 0.55, nsam, n_days),
  nflu_novac = get_counts(200, 50, 0.12, nsam, n_days),
  ve = 0.5,
  lag = 0
)

pop %>%
  mutate(nflu2 = A_to_E + B0_to_F + C_to_F) %>%
  select(nflu2, nflu) %>%
  filter(nflu2 != nflu)

pop %>%
  mutate(averted = nflu_novac - nflu) %>%
  summarise(averted = sum(averted))

plot_distr(200, 50, 0.12, nsam, n_days)



