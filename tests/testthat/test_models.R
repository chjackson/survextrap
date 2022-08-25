postest_test <- function(x, newdata=NULL){
    niter <- 10
    t <- c(1, 5)
    summary(x)
    mean(x, newdata=newdata, niter=niter)
    rmst(x, t=t, newdata=newdata, niter=niter)
    survival(x, newdata=newdata, times=t)
    hazard(x, newdata=newdata, times=t)
    invisible()
}

test_median <- function(mod, vname, value, tol=1e-01){
    expect_equal(summary(mod) %>% filter(variable==vname) %>% pull(median),
                 value, tol=tol)
}

nd <- data.frame(rx = c("Obs", "Lev+5FU"))

test_that("Basic spline model, no covariates",{
    mod <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt")
    test_median(mod, "alpha", -0.578)
    postest_test(mod)
    ## MCMC fit - small sample will give warnings of low ESS etc
    suppressWarnings(modm <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="mcmc",
                                        chains=2, iter=1000, seed=1))
    expect_equal(coef(mod)["alpha"], coef(modm)["alpha"], tol=1e-01)
    ## VB fit. LOO warning
    suppressWarnings(modv <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="vb", loo=FALSE))
    expect_equal(coef(modv)["alpha"], coef(modm)["alpha"], tol=1e-01)
})

test_that("Basic spline model, with covariates",{
    mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
    test_median(mod, "alpha", -0.265)
    postest_test(mod)
    postest_test(mod, nd)
})

test_that("Basic spline model, non-proportional hazards",{
  modnp <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt",
                      nonprop=TRUE, basehaz_ops = list(df=4,degree=2))
  test_median(modnp, "alpha", -0.251)
})


test_that("Weibull model",{
    mod <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt", modelid="weibull")
    test_median(mod, "logscale", 1.6575)
    postest_test(mod)
    mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt", modelid="weibull")
    test_median(mod, "logscale", 1.32)
    postest_test(mod, nd)
    sr <- survival::survreg(Surv(years, status) ~ rx, data=colons, dist="weibull")
    expect_equivalent(coef(sr)["rxLev"], coef(mod)["logtaf_rxLev"], tol=1e-01)
})


test_that("Changing the spline specification",{
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                    basehaz_ops = list(df=4)), "df - degree should be >= 2")
  mod2 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                     basehaz_ops = list(degree=4))
  mod1 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      basehaz_ops = list(degree=1))
  mod0 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      basehaz_ops = list(degree=0))
  expect_error(survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      basehaz_ops = list(degree=-1)), "must be a nonnegative")
  mod0 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      basehaz_ops = list(degree=0, df=2, iknots=1.5))
  expect_equivalent(mod0$basehaz$iknots, 1.5)
})

test_that("Spline prior mean",{
  coef1 <- function(x){summary(x)[summary(x)$variable=="coefs",]$median[1]}
  mod0 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                      basehaz_ops = list(degree=0, df=2, bknots=c(0.1, 4)))
  expect_equivalent(mod0$basehaz$bknots, c(0.1, 4))
  mod01 <-  survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt",
                       basehaz_ops = list(degree=0, df=2, bknots=c(0.1, 4)),
                       coefs_mean = c(0.5, 0.5))
  expect_gt(coef1(mod01), coef1(mod0))
})


test_that("Smoothing standard deviation specifications",{
  mod0 <-  survextrap(Surv(years, status) ~ 1, data=colons, smooth_sd="eb", fit_method="opt",
                      basehaz_ops = list(degree=0, df=2))
  expect_true(mod0$smooth_sd != 1)
  mod1 <-  survextrap(Surv(years, status) ~ 1, data=colons, smooth_sd="bayes", fit_method="opt",
                      basehaz_ops = list(degree=0, df=2))
  expect_equal(mod1$smooth_sd, "bayes")
  mod2 <-  survextrap(Surv(years, status) ~ 1, data=colons, smooth_sd=2, fit_method="opt",
                      basehaz_ops = list(degree=0, df=2))
  expect_equal(mod2$smooth_sd, 2)
})


ndc <- data.frame(x=c(0,1))

test_that("Cure model, no covariates on anything",{
    cmod0 <- survextrap(Surv(t, status) ~ 1, data=curedata, cure=TRUE, fit_method="opt")
    test_median(cmod0, "pcure", 0.574)
    postest_test(cmod0)
})

test_that("Cure model, covariates on noncured model",{
    cmod1 <- survextrap(Surv(t, status) ~ x, data=curedata, cure=TRUE, fit_method="opt")
    test_median(cmod1, "loghr", 0.325)
    postest_test(cmod1, newdata=ndc)
})

test_that("Cure model, covariates on cured fraction",{
    cmod2 <- survextrap(Surv(t, status) ~ 1, data=curedata, cure=~x, fit_method="opt")
    test_median(cmod2, "logor_cure", 0.775)
    postest_test(cmod2, newdata=ndc)
    curedata$x2 <- factor(curedata$x)
    cmod3 <- survextrap(Surv(t, status) ~ 1, data=curedata, cure=~x2, fit_method="opt")
    expect_equal(coef(cmod2)["logor_cure"], coef(cmod3)["logor_cure"], tol=1e-01)
})

test_that("Non-standard model formulae",{
    mod1 <- survextrap(Surv(t, status) ~ 1, data=curedata, cure=~factor(x), fit_method="opt")
    test_median(mod1, "alpha", 2.45)
    survextrap(Surv(t, status) ~ 1, data=curedata, cure=~I(x+1), fit_method="opt")
    survextrap(Surv(t, status) ~ 1, data=curedata, cure=~sqrt(x), fit_method="opt")
    survextrap(Surv(t, status) ~ 1, data=curedata, cure=~splines::bs(x), fit_method="opt")
})

test_that("Relative survival",{
    colonse <- colons
    colonse$bh <- rep(0.01, nrow(colons))
    mod1 <- survextrap(Surv(years, status) ~ 1, data=colonse, backhaz=bh, fit_method="opt")
    colonse$bh <- rep(0.02, nrow(colons))
    mod2 <- survextrap(Surv(years, status) ~ 1, data=colonse, backhaz=bh, fit_method="opt")
    expect_lt(coef(mod2)["alpha"], coef(mod1)["alpha"])

    ext <- data.frame(start=5, stop=10, n=30, r=5,
                      backsurv_start = 0.4, backsurv_stop = 0.3)
    mod1 <- survextrap(Surv(years, status) ~ 1, data=colonse, external=ext, backhaz=bh, fit_method = "opt")

    ext <- data.frame(start=5, stop=10, n=30, r=10,
                      backsurv_start = 0.4, backsurv_stop = 0.3)
    mod2 <- survextrap(Surv(years, status) ~ 1, data=colonse, external=ext, backhaz=bh, fit_method = "opt")
    expect_lt(coef(mod2)["alpha"], coef(mod1)["alpha"])
})
