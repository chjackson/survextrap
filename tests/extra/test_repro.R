test_that("setting the seed", {
  skip_on_cran()
  set.seed(1)
  run1 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                     chains=1, iter=400)
  set.seed(1)
  run2 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                     chains=1, iter=400)
  run2
  # library(testthat)
  expect_equal(summary(run1)$median, summary(run2)$median) # silent = success
  set.seed(2)
  run3 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                     chains=1, iter=400)
  expect_false(all(summary(run1)$median == summary(run3)$median))

  run4 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                     chains=1, iter=400, seed=1)
  run5 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                     chains=1, iter=400, seed=1)
  run6 <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                     chains=1, iter=400, seed=2)
  expect_equal(summary(run4)$median, summary(run5)$median) # success
  expect_false(all(summary(run4)$median == summary(run6)$median))
})
