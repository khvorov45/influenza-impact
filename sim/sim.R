# Simulations like in Tokar 2018
# Arseniy Khvorov
# Created 2019/12/13
# Last edit 2019/12/20

library(tidyverse)
library(lubridate)
library(Rcpp)
library(ggdark) # devtools::install_github("khvorov45/ggdark")

# Directories to be used later
sim_dir <- "sim"

sourceCpp(file.path(sim_dir, "sim.cpp"))
sourceCpp(file.path(sim_dir, "methods.cpp"))

# Functions ===================================================================

# Finds density coefficient to result in a particular overall proportion
find_dens_coef <- function(mean, sd, overall_prop, npop, n_days) {
  overall_prop * npop / sum(dnorm(1:n_days, mean, sd))
}

# Counts with a particular overall proportion
get_counts <- function(mean, sd, overall_prop, npop, n_days) {
  counts <- dnorm(1:n_days, mean, sd) *
    find_dens_coef(mean, sd, overall_prop, npop, n_days)
  as.integer(counts)
}

# Plots a distribution
plot_distr <- function(mean_vac, sd_vac, prop_vac,
                       mean_case, sd_case, prop_case, npop, n_days) {
  tibble(
    day = 1:n_days,
    date = ymd("2010/08/01") + day - 1,
    vacdistr = get_counts(mean_vac, sd_vac, prop_vac, npop, n_days),
    casedistr = get_counts(mean_case, sd_case, prop_case, npop, n_days)
  ) %>%
    pivot_longer(
      contains("distr"), names_to = "distr", values_to = "counts"
    ) %>%
    ggplot(aes(date, counts, lty = distr)) +
    dark_theme_bw(verbose = FALSE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_date(date_breaks = "month") +
    geom_line()
}

# Simulates a population
sim_pop <- function(n_days, init_pop_size, nvac, nflu_novac, ve, lag,
                    seed = sample.int(.Machine$integer.max, 1)) {
  set.seed(seed)
  pop <- sim_pop_cpp(n_days, init_pop_size, nvac, nflu_novac, ve, lag) %>%
    as_tibble() %>%
    mutate(
      date = ymd("2010/08/01") + day - 1,
      averted_true = nflu_novac - nflu
    )
  attr(pop, "seed") <- seed
  attr(pop, "init_pop_size") <- init_pop_size
  attr(pop, "lag") <- lag
  pop
}

month_agg <- function(pop) {
  pop_agg <- pop %>%
    mutate(month = month(date), year = year(date)) %>%
    group_by(year, month) %>%
    summarise(
      vaccinations = sum(A_to_B0 + E_to_F),
      cases = sum(nflu),
      ve = first(ve),
      averted_true = sum(averted_true),
    ) %>%
    ungroup() %>%
    mutate(
      averted_method1 = method1_cpp(cases, vaccinations, ve, nsam) %>%
        pull(averted_method1),
      averted_method2 = method2_cpp(cases, vaccinations, ve, nsam) %>%
        pull(averted_method2),
      averted_method3 = method3_cpp(cases, vaccinations, ve, nsam) %>%
        pull(averted_method3),
      averted_method5 = method5_cpp(cases, vaccinations, ve, nsam) %>%
        pull(averted_method5)
    )
  attr(pop_agg, "seed") <- attr(pop, "seed")
  attr(pop_agg, "init_pop_size") <- attr(pop, "init_pop_size")
  attr(pop_agg, "lag") <- attr(pop, "lag")
  pop_agg
}

# Summarises a population
sum_pop <- function(agg_pop) {
  agg_pop %>%
    summarise(
      averted_true = sum(averted_true),
      averted_method1 = sum(averted_method1),
      averted_method2 = sum(averted_method2),
      averted_method3 = sum(averted_method3),
      averted_method5 = sum(averted_method5),
      init_pop_size = attr(agg_pop, "init_pop_size"),
      lag = attr(agg_pop, "lag"),
      seed = attr(agg_pop, "seed")
    )
}

# Simulate and summarise one
sim_one <- function(n_days, init_pop_size, nvac, nflu_novac, ve, lag,
                    seed = sample.int(.Machine$integer.max, 1)) {
  sourceCpp(file.path(sim_dir, "sim.cpp"))
  sourceCpp(file.path(sim_dir, "methods.cpp"))
  sim_pop(n_days, init_pop_size, nvac, nflu_novac, ve, lag, seed) %>%
    month_agg() %>%
    sum_pop()
}

# Simulate many
sim_many <- function(nsim, n_days, init_pop_size, nvac, nflu_novac, ve, lag,
                     init_seed) {
  map_dfr(
    1:nsim,
    function(ind) sim_one(
      n_days, init_pop_size, nvac, nflu_novac, ve, lag, init_seed + ind - 1
    )
  )
}

# Script ======================================================================

nsim <- 50
nsam <- 1e6
n_days <- 304
#plot_distr(100, 50, 0.55, 190, 35, 0.08, nsam, n_days)
res <- sim_many(
  nsim = nsim,
  n_days = n_days,
  init_pop_size = nsam,
  nvac = get_counts(100, 50, 0.55, nsam, n_days),
  nflu_novac = get_counts(190, 35, 0.12, nsam, n_days),
  ve = 0.48,
  lag = 14,
  init_seed = 1
)

write_csv(
  res, file.path(sim_dir, paste0("res-", nsim, "sims.csv"))
)
