# Simulations like in Tokar 2018
# Arseniy Khvorov
# Created 2019/12/13
# Last edit 2019/01/08

library(tidyverse)
library(lubridate)
library(ggdark) # devtools::install_github("khvorov45/ggdark")
library(impactflu) # devtools::install_github("khvorov45/impactflu")
library(furrr)

plan(multiprocess)

# Directories to be used later
sim_dir <- "sim"

# Functions ===================================================================

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
sim_pop <- function(nsam, vacs, casen, ve, lag, deterministic, start_date,
                    seed = sample.int(.Machine$integer.max, 1)) {
  pop <- sim_ideal(nsam, vacs, casen, ve, lag, deterministic, seed) %>%
    mutate(datestamp = generate_dates(timepoint, start_date, "day"))
  attr(pop, "seed") <- seed
  attr(pop, "lag") <- lag
  attr(pop, "deterministic") <- deterministic
  attr(pop, "nsam") <- nsam
  pop
}

# Monthly aggregate
month_agg <- function(pop) {
  pop_agg <- pop %>%
    mutate(month = month(datestamp), year = year(datestamp)) %>%
    group_by(year, month) %>%
    summarise(
      vaccinations = sum(vaccinations),
      cases = sum(cases),
      ve = mean(ve),
      avert_true = sum(avert),
    ) %>%
    ungroup() %>%
    mutate(
      avert_m1 = method1(nsam, vaccinations, cases, ve) %>% pull(avert),
      avert_m3 = method3(nsam, vaccinations, cases, ve) %>% pull(avert),
    )
  attr(pop_agg, "seed") <- attr(pop, "seed")
  attr(pop_agg, "nsam") <- attr(pop, "nsam")
  attr(pop_agg, "lag") <- attr(pop, "lag")
  attr(pop_agg, "deterministic") <- attr(pop, "deterministic")
  pop_agg
}

# Summarises a monthly aggregate of a population
sum_pop <- function(agg_pop) {
  agg_pop %>%
    summarise(
      avert_true = sum(avert_true),
      avert_m1 = sum(avert_m1),
      avert_m3 = sum(avert_m3),
      nsam = attr(agg_pop, "nsam"),
      lag = attr(agg_pop, "lag"),
      seed = attr(agg_pop, "seed"),
      deterministic = attr(agg_pop, "deterministic")
    )
}

# Simulate and summarise one
sim_one <- function(nsam, vacs, casen, ve, lag, deterministic, start_date,
                    seed = sample.int(.Machine$integer.max, 1)) {
  sim_pop(nsam, vacs, casen, ve, lag, deterministic, start_date, seed) %>%
    month_agg() %>%
    sum_pop()
}

# Simulate many
sim_many <- function(nsim, nsam, vacs, casen, ve, lag, deterministic,
                     start_date, init_seed) {
  future_map_dfr(
    1:nsim,
    function(ind) sim_one(
      nsam, vacs, casen, ve, lag, deterministic, start_date,
      init_seed + ind - 1
    )
  )
}

# Save results
save_res <- function(res, name) {
  write_csv(
    res, file.path(sim_dir, paste0("res-", name, ".csv"))
  )
}

# Script ======================================================================

nsim <- 1e4
nsam <- 1e6
n_days <- 304
start_date = ymd("2017/08/01")

pop <- sim_pop(
  nsam,
  vacs = generate_counts(nsam, n_days, 0.55, mean = 100, sd = 50),
  casen = generate_counts(nsam, n_days, 0.12, mean = 190, sd = 35),
  ve = 0.48,
  lag = 14,
  deterministic = TRUE,
  start_date = start_date
)

pop_month <- month_agg(pop)

pop_summ <- sum_pop(pop_month)

res_deterministic <- sim_one(
  nsam,
  vacs = generate_counts(nsam, n_days, 0.55, mean = 100, sd = 50),
  casen = generate_counts(nsam, n_days, 0.12, mean = 190, sd = 35),
  ve = 0.48,
  lag = 14,
  deterministic = TRUE,
  start_date = start_date
)

save_res(res_deterministic, "deterministic")

res_random <- sim_many(
  nsim, nsam,
  vacs = generate_counts(nsam, n_days, 0.55, mean = 100, sd = 50),
  casen = generate_counts(nsam, n_days, 0.12, mean = 190, sd = 35),
  ve = 0.48,
  lag = 14,
  deterministic = FALSE,
  start_date = start_date,
  init_seed = 1
)

save_res(res_random, paste0(nsim, "sims"))
