# Tests of simulation functions
# Arseniy Khvorov
# Created 2019/12/18
# Last edit 2019/12/18

library(testthat)

test_that("sim_pop works", {
  n_days <- 300
  nsam <- 1e6
  pop <- sim_pop(
    n_days = n_days,
    init_pop_size = nsam,
    nvac = get_counts(100, 50, 0.55, nsam, n_days),
    nflu_novac = get_counts(200, 50, 0.12, nsam, n_days),
    ve = 0.5,
    lag = 1
  )
  pop_wrong <- pop %>%
    mutate(nflu2 = A_to_E + B0_to_F + C_to_F + B1_to_F) %>%
    select(nflu2, nflu) %>%
    filter(nflu2 != nflu)
  expect_equal(nrow(pop_wrong), 0)
})

