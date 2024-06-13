extdat <- data.frame(start = c(5, 10, 15, 20),
                     stop =  c(10, 15, 20, 25),
                     n = c(100, 100, 100, 100),
                     r = c(50, 40, 30, 20))

test_that("External data gives precise extrapolation",{
  nde_mod <- survextrap(Surv(years, status) ~ 1, data=colons,
                        chains=1, external = extdat,
                        add_knots=c(4, 10, 25), fit_method="opt")
  expect_true(hazard(nde_mod, t=20)$upper < 0.5)
  expect_true(all(c(4,10,25) %in% nde_mod$mspline$knots))
})

test_that("External data with covariates",{
  extdat$sex = 1
  nde_mod1 <- survextrap(Surv(years, status) ~ sex, data=colons,
                         chains=1, external = extdat,
                         add_knots=c(4, 10, 25), fit_method="opt")
  haz1 <- hazard(nde_mod1, newdata=list(sex=c(0,1)), t=20)
  extdat$sex = 0
  nde_mod0 <- survextrap(Surv(years, status) ~ sex, data=colons,
                         chains=1, external = extdat,
                         add_knots=c(4, 10, 25), fit_method="opt")
  haz0 <- hazard(nde_mod0, newdata=list(sex=c(0,1)), t=20)
  ## more precise estimate given external data with the desired covariate value
  expect_true(haz1$upper[haz1$sex==1] < haz0$upper[haz0$sex==1])
})

test_that("External data with nonproportional hazards",{
  extdat$sex = 1
  expect_no_error({
    nde_mod1 <- survextrap(Surv(years, status) ~ sex, nonprop=TRUE, data=colons,
                           chains=1, external = extdat,
                           add_knots=c(4, 10, 25), fit_method="opt")
    nde_mod1
  })
})

test_that("External data with cure and covariates",{
  extdat$sex = 1
  expect_no_error({
    nde_mod1 <- survextrap(Surv(years, status) ~ 1, cure=~sex, data=colons,
                           chains=1, external = extdat,
                           add_knots=c(4, 10, 25), fit_method="opt")
    nde_mod1
  })
})

test_that("external data MCMC fit",{
  skip_on_cran()
  extdat$sex = 1
  suppressWarnings({
    nde_mod1 <- survextrap(Surv(years, status) ~ sex, data=colons,
                           chains=1, external = extdat,
                           add_knots=c(4, 10, 25), fit_method="mcmc", iter=1000)
    expect_true(is.numeric(nde_mod1$loo$estimates["looic","Estimate"]))
  })
})

test_that("No individual-level data",{
  nde_mod1 <- survextrap(~1, external = extdat,
                         mspline = list(df=5), fit_method="opt")
  nde_mod1$mspline
  plot_hazsurv(nde_mod1)
  expect_warning(nde_mod1 <- survextrap(Surv(time, event) ~1, external = extdat,
                                        mspline = list(df=5), fit_method="opt"),
                 "`formula` has a left-hand side, but `data` not supplied")
  expect_error(survextrap(external = extdat), "argument \"formula\" is missing")
})

extdatc <- rbind(cbind(extdat, treat="no"),
                 cbind(extdat, treat="yes"))
extdatc$treat <- factor(extdatc$treat)

test_that("No individual-level data, covariates",{
  nde_mod2 <- survextrap(~treat, external = extdatc,
                         mspline = list(df=5), fit_method="opt")
  expect_false(nde_mod2$indiv)
})
