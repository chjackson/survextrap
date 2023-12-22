alpha <- -0.35
coefs <- c(0.0012, 0.029, 0.1, 0.15, 0.23, 0.0074, 0.052, 0.1)
knots <- c(0.36, 0.62, 0.91, 1.2, 1.54, 2.02, 3)

test_that("wane functions",{
  pnw <- psurvmspline(q=c(10,14,18), alpha=alpha, coefs=coefs, knots=knots)
  pw <- psurvmspline_wane(q=10, alpha1=alpha, alpha0=-0.2, coefs1=coefs, coefs0=coefs,
                          knots=knots, wane_period=c(11, 15))
  hnw <- hsurvmspline(x=c(10,14,18), alpha=alpha, coefs=coefs, knots=knots)
  expect_equal(hnw[1],hnw[3])
  hw <- hsurvmspline_wane(x=c(10,14,18), alpha1=alpha, alpha0=-0.2, coefs1=coefs, coefs0=coefs,
                          knots=knots, wane_period=c(11, 15))
  dnw <- dsurvmspline(x=10, alpha=alpha, coefs=coefs, knots=knots)
  dw <- dsurvmspline_wane(x=10, alpha1=alpha, alpha0=-0.2, coefs1=coefs, coefs0=coefs,
                          knots=knots, wane_period=c(11, 15))
  expect_equal(pnw[1], pw[1])
  expect_equal(dnw[1], dw[1])
  expect_equal(hnw[1], hw[1])
  expect_lt(hnw[2], hw[2])
  expect_lt(hnw[3], hw[3])
})

test_that("wane functions: vectorisation over parameters",{
  h1 <- Hsurvmspline_wane(x=c(8, 10), alpha1=-1, alpha0=-2,
                          coefs1=coefs, coefs0=coefs, knots=knots, wane_period=c(9, 11))
  h2 <- Hsurvmspline_wane(x=c(8, 10), alpha1=-3, alpha0=-4,
                          coefs1=coefs, coefs0=coefs, knots=knots, wane_period=c(9, 11))
  h12 <- Hsurvmspline_wane(x=c(8, 10), alpha1=c(-1,-3), alpha0=c(-2,-4),
                           coefs1=coefs, coefs0=coefs, knots=knots, wane_period=c(9, 11))
  expect_equal(h12, c(h1[1], h2[2]))
})

test_that("wane functions: vectorisation errors",{
  expect_error(Hsurvmspline_wane(x=c(0, 8, 10, 12), alpha1=rep(alpha,2), alpha0=alpha*1.01,
                               coefs1=coefs, coefs0=coefs, knots=knots, wane_period=c(9, 11)),
             "lengths of `alpha0`")
})

mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")

test_that("default newdata0",{
  nd <- data.frame(rx = c("Obs", "Lev+5FU"))
  nd0 <- data.frame(rx = c("Obs", "Obs"))
  rd <- rmst(mod, wane_period=c(2,3), t=7, niter=5, wane_nt=2)
  rn <- rmst(mod, newdata0=nd0, newdata=nd, wane_period=c(2,3), t=7, niter=5, wane_nt=2)
  expect_equal(rd$median[c(1,3)], rn$median, tol=1e-06)
})

test_that("rmst wane",{
  rd <- rmst(mod, t=2, niter=5)
  rd2 <- rmst(mod, t=2, wane_period=c(5,7), niter=5)
  expect_equal(rd$median, rd2$median)
  rd <- rmst(mod, t=5, niter=5)
  rd2 <- rmst(mod, t=5, wane_period=c(3,5), niter=5)
  expect_lt(rd2$median[2], rd$median[2])
})

