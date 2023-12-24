test_that("loo",{
  skip_on_cran()
  set.seed(1)
  suppressWarnings({
    mod <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                      chains=1, iter=1000)
    expect_equal(mod$loo$estimates["looic","Estimate"], 426.7, tol=1e-02)
    modr <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="mcmc",
                       chains=1, iter=1000)
    expect_equal(modr$loo$estimates["looic","Estimate"], 425.01, tol=1e-02)
    expect_lt(modr$loo$estimates["looic","Estimate"],
              mod$loo$estimates["looic","Estimate"])
  })
})

test_that("loo with pcure",{
  skip_on_cran()
  set.seed(1)
  suppressWarnings({
    cmod0 <- survextrap(Surv(t, status) ~ 1, data=curedata, cure=TRUE,
                        fit_method="mcmc", chains=1, iter=1000)
    expect_equal(cmod0$loo$estimates["looic","Estimate"], 464.6158, tol=1e-02)
  })
})

test_that("loo with no events / censored",{
  skip_on_cran()
  set.seed(1)
  cnoev <- colons[colons$status==0,]
  cnocens <- colons[colons$status==1,]
  suppressWarnings({
    mod <- survextrap(Surv(years, status) ~ 1, data=cnoev, fit_method="mcmc",
                      chains=1, iter=1000)
    expect_equal(mod$loo$estimates["looic","Estimate"], 0.0647417, tol=1e-02)
    mod <- survextrap(Surv(years, status) ~ 1, data=cnocens, fit_method="mcmc",
                      chains=1, iter=1000)
    expect_equal(mod$loo$estimates["looic","Estimate"], 169.6696, tol=1e-02)
  })
})

test_that("loo with external",{
  skip_on_cran()
  extdat <- data.frame(start = c(5, 10, 15, 20),
                       stop =  c(10, 15, 20, 25),
                       n = c(100, 100, 100, 100),
                       r = c(50, 40, 30, 20))
  set.seed(1)
  suppressWarnings({
    nde_mod <- survextrap(Surv(years, status) ~ 1, data=colons,
                          external = extdat, add_knots=c(4, 10, 25),
                          fit_method="mcmc", chains=1, iter=1000)
    expect_equal(nde_mod$loo$estimates["looic","Estimate"], 426.8256, tol=1e-03)

    extdat$sex = 1
    nde_mod1 <- survextrap(Surv(years, status) ~ sex, data=colons,
                           external = extdat, add_knots=c(4, 10, 25),
                           fit_method="mcmc", chains=1, iter=1000)
    expect_equal(nde_mod1$loo$estimates["looic","Estimate"], 428.6499, tol=1e-03)
  })
})
