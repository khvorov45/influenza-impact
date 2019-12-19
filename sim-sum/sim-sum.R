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
    contains("averted_method"), names_to = "method", values_to = "averted"
  ) %>%
  mutate(
    method = str_replace(method, "averted_method", ""),
    bias = averted - averted_true,
    bias_prop = bias / averted_true
  ) %>%
  group_by(method, init_pop_size, lag) %>%
  summarise(
    bias_mean = mean(bias),
    bias_prop_mean = mean(bias_prop),
    averted_sd = sd(averted),
  ) %>%
  ungroup()

write_csv(summ, file.path(sim_sum_dir, "sum-10000sims.csv"))
