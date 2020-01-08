# Summary of sim
# Arseniy Khvorov
# Created 2019/12/19
# Last edit 2019/12/19

library(tidyverse)

# Directories to be used later
sim_dir <- "sim"
sim_sum_dir <- "sim-sum"

# Script ======================================================================

res <- read_csv(file.path(sim_dir, "res-10000sims.csv"))

summ <- res %>%
  pivot_longer(
    contains("avert_m"), names_to = "method", values_to = "avert"
  ) %>%
  mutate(
    method = str_replace(method, "avert_m", ""),
    bias = avert - avert_true,
    bias_prop = bias / avert_true
  ) %>%
  group_by(method, nsam, lag) %>%
  summarise(
    bias_mean = mean(bias),
    bias_prop_mean = mean(bias_prop),
    averted_sd = sd(avert),
  ) %>%
  ungroup()

write_csv(summ, file.path(sim_sum_dir, "sum-10000sims.csv"))
