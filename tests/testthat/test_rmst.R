## Integration methods for RMST and related functions

mod <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt")

test_that("RMST integration methods at vector of times",{
  intnew <- rmst(mod, t=c(1,7), niter=5, method="gl")
  intold <- rmst(mod, t=c(1,7), niter=5, method="adaptive")
  expect_equal(intnew, intold, tolerance=1e-04)
})

modc <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")

test_that("RMST integration methods with covariates",{
  intnew <- rmst(modc, t=c(1,7), niter=5, method="gl")
  intold <- rmst(modc, t=c(1,7), niter=5, method="adaptive")
  expect_equal(intnew, intold, tolerance=1e-04)
})

test_that("RMST integration methods with waning",{
  intnew <- rmst(modc, t=c(1,7), niter=5, wane_period = c(6, 9), method="gl")
  intold <- rmst(modc, t=c(1,7), niter=5, wane_period = c(6, 9), method="adaptive")
  expect_equal(intnew, intold, tolerance=1e-04)
})

test_that("RMST integration methods are different but close enough",{
  intold <- rmst(mod, t=7, niter=5, method="adaptive")
  intnew <- rmst(mod, t=7, niter=5, method="gl")
  intnew10 <- rmst(mod, t=7, niter=5, method="gl", gl_nodes=10)
  expect_false(isTRUE(all.equal(intold$median, intnew$median, tolerance=1e-10)))
  expect_false(isTRUE(all.equal(intnew$median, intnew10$median, tolerance=1e-10)))
  expect_equal(intold$median, intnew$median, tolerance=1e-05)
  intnew <- rmst_generic(plnorm, t=100, meanlog=0.1, sdlog=0.2, method="gl")
  intold <- rmst_generic(plnorm, t=100, meanlog=0.1, sdlog=0.2, method="adaptive")
  expect_equal(intnew, intold, tolerance=1e-04)
})

test_that("rmst_generic matches analytic RMST", {
  rate <- 0.1
  ref <- 1/rate * pexp(5, rate=rate)
  intold <- rmst_generic(pexp, t=5, rate=rate, method="adaptive")
  intnew <- rmst_generic(pexp, t=5, rate=rate, method="gl")
  expect_equal(intnew, ref)
  expect_equal(intold, ref)
})

test_that("rmst_generic works with different argument formats", {
  intnew <- rmst_generic(pexp, t=5, rate=c(0.3, 0.4), method="gl")
  ref <- 1/c(0.3, 0.4) * pexp(5, rate=c(0.3, 0.4))
  expect_equal(intnew, ref)
  intnew <- rmst_generic(pexp, t=c(5,10), rate=0.4, method="gl")
  ref <- 1/0.4 * pexp(c(5,10), rate=0.4)
  expect_equal(intnew, ref)
  ref <- (1/0.4 * (pexp(c(5,10), rate=0.4) - pexp(c(0, 1), rate=0.4))) / (1 - pexp(c(0,1),rate=0.4))
  intnew <- rmst_generic(pexp, start=c(0, 1), t=c(5,10), rate=0.4, method="gl")
  expect_equal(intnew, ref)
})

test_that("errors in RMST specification",{
  expect_error(rmst(mod, t=7, niter=5, method="wibble"), "unknown integration `method`")
})

test_that("unrestricted mean survival time uses adaptive method",{
  intold <- mean(mod, niter=5, method="adaptive")
  intnew <- mean(mod, niter=5, method="gl")
  expect_equal(intold$median, intnew$median)
})


modbh <- survextrap(Surv(years, d) ~ 1, data = cetux,
                    backhaz=cetux_bh, fit_method = "opt")
modbhc <- survextrap(Surv(years, d) ~ treat, data = cetux,
                    backhaz=cetux_bh, fit_method = "opt")

test_that("RMST integration methods with backhaz",{
  intnew <- rmst(modbh, t=c(1,7), niter=1, method="gl")
  intold <- rmst(modbh, t=c(1,7), niter=1, method="adaptive")
  expect_equal(intold$median, intnew$median, tol=1e-05)
})
