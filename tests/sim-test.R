# Tests of simulation functions
# Arseniy Khvorov
# Created 2019/12/18
# Last edit 2019/12/23

library(testthat)

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
pop_sum <- month_agg(pop)

test_that("sim_pop works", {
  pop_check <- pop %>%
    mutate(
      pflu2 = nflu_novac / lag(uninf_novac, default = 1),
      nflu2 = A_to_E + B0_to_F + C_to_F + B1_to_F,
      uninf_novac2 = lag(uninf_novac, default = nsam) - nflu_novac,
      avert2 = nflu_novac - nflu,
      pvac2 = nvac / (lag(A, default = nsam) + lag(E, default = 0)),
      A2 = lag(A, default = nsam) - A_to_B0 - A_to_E,
      B1_2 = lag(B0, default = 0) - B0_to_F,
      C2 = lag(C, default = 0) - C_to_F + B1_to_C,
      D2 = lag(D, default = 0) + B1_to_D,
      E2 = lag(E, default = 0) + A_to_E - E_to_F,
      F2 = lag(`F`, default = 0) + B0_to_F + B1_to_F + C_to_F + E_to_F
    ) %>%
    select(
      pflu2, pflu, nflu2, nflu, uninf_novac2, uninf_novac, avert2, averted_true,
      pvac2, pvac, A2, A, B1_2, B1, C2, C, D2, D, E2, E, F2, `F`
    )
  with(pop_check, {
    expect_equal(pflu2, pflu)
    expect_equal(nflu2, nflu)
    expect_equal(uninf_novac2, uninf_novac)
    expect_equal(avert2, averted_true)
    expect_equal(pvac2, pvac)
    expect_equal(A2, A)
    expect_equal(B1_2, B1)
    expect_equal(C2, C)
    expect_equal(D2, D)
    expect_equal(E2, E)
    expect_equal(F2, `F`)
  })
})

test_that("method1 works", {
  m1 <- method1_cpp(pop_sum$cases, pop_sum$vaccinations, pop_sum$ve, nsam)
  m1_check <- as_tibble(m1) %>%
    mutate(
      pvac = vaccinations / nsam,
      vc_lag2 = (pvac + lag(pvac, default = 0)) / 2,
      pops2 = (lag(pops, default = nsam) - lag(cases, default = 0)) *
        (1 - vc_lag * ve),
      pops2 = as.integer(pops2),
      pflu2 = cases / pops,
      popn2 = lag(popn, default = nsam - casen[[1]]) - lag(casen, default = 0),
      popn2 = as.integer(popn2),
      casen2 = pflu * popn,
      casen2 = floor(casen2) %>% as.integer(),
      averted2 = casen - cases
    ) %>%
    select(
      vc_lag, vc_lag2, pops2, pops, pflu2, pflu, popn2, popn, casen2, casen,
      averted2, averted_method1
    )
  with(m1_check, {
    expect_equal(vc_lag2, vc_lag)
    expect_equal(pops2, pops)
    expect_equal(pflu2, pflu)
    expect_equal(popn2, popn)
    expect_equal(casen2, casen)
    expect_equal(averted2, averted_method1)
  })
})

test_that("method2 works", {
  m2 <- method2_cpp(pop_sum$cases, pop_sum$vaccinations, pop_sum$ve, nsam)
  m2_check <- as_tibble(m2) %>%
    mutate(
      pvac = vaccinations / nsam,
      vef2 = lag(pop, default = nsam) * pvac * ve,
      vef2 = as.integer(vef2),
      pop2 = lag(pop, default = as.integer(nsam)) - cases,
      pops2 = lag(pops, default = nsam) - cases - vef,
      pflu2 = cases / lag(pops, default = nsam),
      popn2 = lag(popn, default = nsam) - lag(casen, default = 0),
      popn2 = as.integer(popn2),
      casen2 = pflu * popn,
      casen2 = floor(casen2) %>% as.integer(),
      averted2 = casen - cases
    ) %>%
    select(
      vef2, vef, pop2, pop, pops2, pops, pflu2, pflu, popn2, popn, casen2,
      casen, averted2, averted_method2
    )
  with(m2_check, {
    expect_equal(vef2, vef)
    expect_equal(pop2, pop)
    expect_equal(pops2, pops)
    expect_equal(pflu2, pflu)
    expect_equal(popn2, popn)
    expect_equal(casen2, casen)
    expect_equal(averted2, averted_method2)
  })
})
