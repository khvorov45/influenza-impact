# Simulations like in Tokar 2018
# Arseniy Khvorov
# Created 2019/12/13
# Last edit 2019/12/13

library(tidyverse)
library(Rcpp)
library(extraDistr)

# Directories to be used later
sim_dir <- "sim"

sourceCpp(file.path(sim_dir, "sim.cpp"))

# Functions ===================================================================

# Simulate a population on a particular day
sim_pop <- function(pop_prev, pvac, pflu, ve, vac_lag_ago = NULL) {
  pop <- pop_prev %>%
    mutate(
      vac = if_else(vac == 0 & inf == 0, rbern(n(), pvac), vac),
      inf = if_else(sus == 1, rbern(n(), pflu), inf),
      sus = case_when(
        #sus == 0 ~ 0,
        inf == 1 ~ 0,
        id %in% vac_lag_ago ~ rbern(n(), 1 - ve),
        TRUE ~ sus
      )
    )
  attr(pop, "pvac") <- pvac
  attr(pop, "pflu") <- pflu
  attr(pop, "ve") <- ve
  pop
}

# Summarise a particular day's population
sum_pop <- function(pop, day) {
  pop %>%
    summarise(
      day = day,
      unvac_uninf = sum(vac == 0 & inf == 0),
      unvac_inf = sum(vac == 0 & inf == 1),
      vac_uninf_sus = sum(vac == 1 & inf == 0 & sus == 1),
      vac_uninf_imm = sum(vac == 1 & inf == 0 & sus == 0),
      vac_inf = sum(vac == 1 & inf == 1),
      pvac = attr(pop, "pvac"),
      pflu = attr(pop, "pflu"),
      ve = attr(pop, "ve")
    )
}

# Simulate a population after a sequence of days
sim_pop_final <- function(nsam, ndays, pvac, pflu, ve, lag) {
  pop_init <- tibble(
    id = 1:nsam,
    vac = rep(0, nsam),
    inf = rep(0, nsam),
    sus = rep(1, nsam)
  )
  attr(pop_init, "pvac") <- 0
  attr(pop_init, "pflu") <- 0
  attr(pop_init, "ve") <- ve
  res <- list(sum_pop(pop_init, 1))
  pops <- list(pop_init)
  for (i in 2:ndays) {
    if (i > lag + 1) {
      vac_lagm1_ago <- pops[[i - lag - 1]] %>% filter(vac == 1) %>% pull(id)
      vac_lag_ago <- pops[[i - lag]] %>% filter(vac == 1) %>% pull(id)
      vac_lag_ago <- vac_lag_ago[!vac_lag_ago %in% vac_lagm1_ago]
      still_uninf <- pops[[i - 1]] %>% filter(inf == 0) %>% pull(id)
      vac_lag_ago <- intersect(vac_lag_ago, still_uninf)
    } else vac_lag_ago <- NULL
    pops[[i]] <- sim_pop(
      pops[[i - 1]], pvac[[i - 1]], pflu[[i - 1]], ve, vac_lag_ago
    )
    res[[i]] <- sum_pop(pops[[i]], i)
  }
  bind_rows(res)
}

# Script ======================================================================

nsam <- 1e5
ndays <- 200
res <- sim_pop_final(
  nsam, ndays, c(rep(0.05, ndays - 2), 0), rep(0.05, ndays - 1), 0.4, 1
)





