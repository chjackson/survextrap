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

test_that("Weibull model",{
    mod <- survextrap(Surv(years, status) ~ 1, data=colons, fit_method="opt", modelid="weibull")
    test_median(mod, "alpha", 1.6575)
    postest_test(mod)
    mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt", modelid="weibull")
    test_median(mod, "alpha", 1.32)
    postest_test(mod, nd)
    sr <- survival::survreg(Surv(years, status) ~ rx, data=colons, dist="weibull")
    expect_equivalent(coef(sr)["rxLev"], coef(mod)["loghr_rxLev"], tol=1e-01)
})


### Changing various settings in basic model
### Internal knots, boundary knots
### Spline prior mean
### EB for smooth variance
### Fixed smooth variance
### With external data



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


## TODO test error handling for if wrong covariate names supplied in newdata
## And a lot more errors

## Specific tests for post estimation functions
