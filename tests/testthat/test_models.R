postest_test <- function(x, newdata=NULL){
    niter <- 10
    t <- c(1, 5)
    summary(x)
    rmst(x, t=t, newdata=newdata, niter=niter)
    survival(x, newdata=newdata, t=t)
    hazard(x, newdata=newdata, t=t)
    invisible()
}

test_median <- function(mod, vname, value, tol=1e-01){
    expect_equal(summary(mod) %>% filter(variable==vname) %>% pull(median) %>% as.numeric(),
                 value, tol=tol)
}

nd <- data.frame(rx = c("Obs", "Lev+5FU"))

test_that("Basic spline model, no covariates",{
    mod <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt")
    test_median(mod, "alpha", -0.578)
    postest_test(mod)
    ## MCMC fit - small sample will give warnings of low ESS etc
    suppressWarnings(modm <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                                        chains=1, iter=400, seed=1))
    expect_equal(coef(mod)["alpha"], coef(modm)["alpha"], tol=1e-01)
    expect_true(is.numeric(mod$prior_sample$sample(nsim=4)$alpha))
    expect_true(is.numeric(mod$prior_sample$haz_const()["50%","haz"]))
})

test_that("Basic spline model, with covariates",{
    mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
    test_median(mod, "alpha", -0.265)
    postest_test(mod)
    postest_test(mod, nd)
    mod$prior_sample$sample(nsim=4)
    mod$prior_sample$sample(nsim=4, newdata=data.frame(rx="Lev"))
})

test_that("Basic spline model, non-proportional hazards",{
  modnp <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt",
                      nonprop=TRUE, mspline = list(df=4,degree=2,bsmooth=FALSE))
  test_median(modnp, "alpha", -0.251)
  expect_equal(survival(modnp, newdata=nd, t=2)$median, c(0.5, 0.7), tol=1e-01)
  expect_equal(hazard(modnp, newdata=nd, t=2)$median, c(0.2, 0.1), tol=1e-01)
  modnp <- survextrap(Surv(years, status) ~ rx + sex + obstruct,
                      data=colons, fit_method="opt",
                      nonprop=~sex + obstruct,
                      mspline = list(df=4,degree=2,bsmooth=FALSE))
  expect_equal(summary(modnp) %>% filter(variable=="hrsd") %>% pull(term),
               c("sex", "obstruct"))
})


test_that("Changing the spline specification",{
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                          mspline = list(df=4, bsmooth=FALSE)),
               "df - degree should be >= 2")
  expect_no_error(survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                          mspline = list(df=4, bsmooth=TRUE)))
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                          mspline = list(df=2, bsmooth=TRUE)),
               "df - degree should be >= 0")
  mod2 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                     mspline = list(degree=4, bsmooth=FALSE))
  mod1 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      mspline = list(degree=1, bsmooth=FALSE))
  mod0 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      mspline = list(degree=0, bsmooth=FALSE))
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      mspline = list(degree=-1, bsmooth=FALSE)),
               "must be a nonnegative")
  mod0 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      mspline = list(degree=0, knots=c(1.5, 3),
                                     bsmooth=FALSE))
  expect_equivalent(mod0$mspline$knots, c(1.5, 3))

  mspline <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=6, add_knots=20)
  expect_equal(mspline$knots[[2]], 1.04)

})

test_that("Random walk priors",{
  expect_no_error({
    modr <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                       smooth_model = "random_walk")
    rxnphr_mod <- survextrap(Surv(years, status) ~ rx, data=colons,
                            nonprop=TRUE, fit_method = "opt",
                            smooth_model="random_walk")
  })
})

test_that("Spline prior mean",{
  coef1 <- function(x){summary(x)[summary(x)$variable=="coefs",]$median[1]}
  mod0 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      mspline = list(degree=1, knots=c(2,4), bsmooth=FALSE))
  expect_equivalent(mod0$mspline$knots, c(2, 4))
  mod01 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                       mspline = list(degree=1, knots=c(2, 4), bsmooth=FALSE),
                       coefs_mean = c(0.98, 0.01, 0.01))
  expect_gt(coef1(mod01), coef1(mod0))
})

test_that("Smoothing standard deviation specifications",{
  mod0 <-  survextrap(Surv(years, status) ~ 1, data=colons, hsd="eb", fit_method="opt",
                      mspline = list(degree=0, df=2, bsmooth=FALSE))
  expect_true(mod0$hsd != 1)
  mod1 <-  survextrap(Surv(years, status) ~ 1, data=colons, hsd="bayes", fit_method="opt",
                      mspline = list(degree=0, df=2, bsmooth=FALSE))
  expect_equal(mod1$hsd, "bayes")
  mod2 <-  survextrap(Surv(years, status) ~ 1, data=colons, hsd=2, fit_method="opt",
                      mspline = list(degree=0, df=2, bsmooth=FALSE))
  expect_equal(mod2$hsd, 2)
})


ndc <- data.frame(x=c(0,1))

test_that("Cure model, no covariates on anything",{
    cmod0 <- survextrap(Surv(t, status) ~ 1, mspline=list(bsmooth=FALSE),
                        data=curedata, cure=TRUE, fit_method="opt")
    test_median(cmod0, "pcure", 0.574)
    postest_test(cmod0)
})

test_that("Cure model, covariates on noncured model",{
  set.seed(1)
  cmod1 <- survextrap(Surv(t, status) ~ x, data=curedata, cure=TRUE,
                      mspline=list(bsmooth=FALSE), fit_method="opt")
  test_median(cmod1, "loghr", 0.3)
  postest_test(cmod1, newdata=ndc)
})

test_that("Cure model, covariates on cured fraction",{
    cmod2 <- survextrap(Surv(t, status) ~ 1, mspline=list(bsmooth=FALSE),
                        data=curedata, cure=~x, fit_method="opt")
    test_median(cmod2, "logor_cure", 0.775)
    postest_test(cmod2, newdata=ndc)
    curedata$x2 <- factor(curedata$x)
    cmod3 <- survextrap(Surv(t, status) ~ 1, data=curedata, cure=~x2, fit_method="opt")
    expect_equal(coef(cmod2)["logor_cure"], coef(cmod3)["logor_cure"], tol=1e-01)
})

test_that("Non-standard model formulae",{
  mod1 <- survextrap(Surv(t, status) ~ 1, data=curedata, cure=~factor(x),
                     mspline=list(bsmooth=FALSE), fit_method="opt")
  haz <- hazard(mod1, t=3, niter=1)
  expect_equal(haz$median[haz$x==0], 0.0895, tol=0.1)
  test_median(mod1, "logor_cure", 0.775)
  curedata$xf <- factor(curedata$x)
  mod1 <- survextrap(Surv(t, status) ~ 1, data=curedata, cure=~xf, fit_method="opt")
  survextrap(Surv(t, status) ~ 1, data=curedata, cure=~I(x+1), fit_method="opt")
  survextrap(Surv(t, status) ~ 1, data=curedata, cure=~sqrt(x), fit_method="opt")
  survextrap(Surv(t, status) ~ 1, data=curedata, cure=~splines::bs(x), fit_method="opt")
})

test_that("Relative survival models specified through a variable in the data",{
    colonse <- colons
    colonse$bh <- rep(0.01, nrow(colons))
    mod1 <- survextrap(Surv(years, status) ~ 1, data=colonse, backhaz="bh", fit_method="opt")
    colonse$bh <- rep(0.02, nrow(colons))
    mod2 <- survextrap(Surv(years, status) ~ 1, data=colonse, backhaz="bh", fit_method="opt")
    expect_lt(coef(mod2)["alpha"], coef(mod1)["alpha"])

    ext <- data.frame(start=5, stop=10, n=30, r=5,
                      backsurv_start = 0.4, backsurv_stop = 0.3)
    mod1 <- survextrap(Surv(years, status) ~ 1, data=colonse, external=ext, backhaz="bh", fit_method = "opt")

    ext <- data.frame(start=5, stop=10, n=30, r=10,
                      backsurv_start = 0.4, backsurv_stop = 0.3)
    mod2 <- survextrap(Surv(years, status) ~ 1, data=colonse, external=ext, backhaz="bh", fit_method = "opt")
    expect_lt(coef(mod2)["alpha"], coef(mod1)["alpha"])
})

test_that("Relative survival models specified through a background hazard data frame",{
  bh <- data.frame(hazard = c(0.01, 0.02, 0.03), time=c(0, 5, 10))
  mod1 <- survextrap(Surv(years, status) ~ 1, data=colons, backhaz=bh, fit_method="opt")
  rmst(mod1, t=10, niter=10)
  plot_hazard(mod1, niter=10, tmax=20)
  ## extrapolation is accounting for increase in background at times 5 and 10
  ## but uncertainty also because cause-specific hazard is being extrapolated
  haz <- hazard(mod1, t=c(5,10))
  expect_lt(haz$median[1], haz$median[2])
})

test_that("Cure model coupled with a background hazard data frame",{
  expect_no_error({
    cmod0 <- survextrap(Surv(t, status) ~ 1, data=curedata, cure=TRUE, fit_method="opt")
    plot_hazard(cmod0, niter=20)

    ## Use cure model for short term, and background for long term
    bh <- data.frame(hazard = c(0.01, 0.05, 0.1, 0.5), time=c(0, 5, 7, 10))
    cmod1 <- survextrap(Surv(t, status) ~ 1, data=curedata, cure=TRUE,
                        backhaz=bh, fit_method="opt")
    plot_hazard(cmod1, tmax=12, niter=50) + coord_cartesian(ylim=c(0,1))
    plot_survival(cmod1, tmax=20, niter=50)
  })
})

test_that("Cure and relative survival with MCMC",{
  skip_on_cran()
  suppressWarnings({
    cmod0 <- survextrap(Surv(t, status) ~ 1, mspline=list(bsmooth=FALSE),
                        data=curedata, cure=TRUE, fit_method="mcmc",chains=1, iter=1000)
    expect_true(is.numeric(cmod0$loo$estimates["looic","Estimate"]))
    colonse <- colons
    colonse$bh <- rep(0.01, nrow(colons))
    mod1 <- survextrap(Surv(years, status) ~ 1, data=colonse, backhaz="bh",
                       fit_method="mcmc", chains=1, iter=1000)
    expect_true(is.numeric(mod1$loo$estimates["looic","Estimate"]))
  })
})
