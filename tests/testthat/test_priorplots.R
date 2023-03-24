
test_that("spline demo plots", {
  expect_error(mspline_plotdata(iknots = c(1,2,3), bknots = c(0, 5), df=5, coefs=c(1, rep(10, 3), 1)),
               "length of `coefs` is 5, should be 7")
})

sp <- list(iknots=1:3, bknots=c(0,5), degree=2)

test_that("prior hazard parameter samples", {
  expect_error({ # ?? any sensible expectations here?
    ## no covariates
    prior_sample(mspline = sp, 
                 coefs_mean = NULL,
                 prior_smooth = p_gamma(2,1),
                 prior_loghaz = p_normal(0,1),
                 nsim = 100)

    ## prop haz model, 1 cov
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_smooth = p_gamma(2,1),
                 prior_loghaz = p_normal(0,1),
                 x = list(ncovs=1, xnames="treatment"),
                 X = list(treatment=1),
                 nsim = 100)
    ## prop haz model, 2 covs
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_smooth = p_gamma(2,1),
                 prior_loghaz = p_normal(0,1),
                 x = list(ncovs=2, xnames=c("treatment","age")),
                 X = list(treatment=1, age=60),
                 nsim = 100)

    ## nonprop haz model, 1 cov
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_smooth = p_gamma(2,1),
                 prior_loghaz = p_normal(0,1),
                 nonprop = TRUE,
                 x = list(ncovs=1, xnames="treatment"),
                 X = list(treatment=1),
                 nsim = 100)
    ## nonprop haz model, 2 covs
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_smooth = p_gamma(2,1),
                 prior_loghaz = p_normal(0,1),
                 nonprop = TRUE,
                 x = list(ncovs=1, xnames="treatment"),
                 X = list(treatment=1),
                 nsim = 100)
  },NA)
})


test_that("prior hazard function samples", {
  expect_error({
    pdf <- mspline_priorpred(iknots=1:3, bknots=c(0,5), degree=2,
                                prior_loghaz = p_normal(0,1),
                                prior_smooth = p_gamma(10, 10),
                                tmin=0, tmax=10, nsim=10)
    ggplot(pdf, aes(x=time, y=haz, group=rep)) + geom_line() + ylim(0,2)

    ## prop haz model
    pdf <- mspline_priorpred(iknots=1:3, bknots=c(0,5), degree=2,
                                prior_loghaz = p_normal(0,1),
                                prior_smooth = p_gamma(10, 10),
                                x = list(ncovs=1, xnames="treatment"),
                                X = list(treatment=-5),
                                tmin=0, tmax=10, nsim=10)
    ggplot(pdf, aes(x=time, y=haz, group=rep)) + geom_line() + ylim(0,2)

    ## nonprop haz model
    p <- plot_mspline_priorpred(iknots=1:3, bknots=c(0,5), df=5, degree=2,
                                prior_loghaz = p_normal(0,1),
                                prior_smooth = p_gamma(10, 10),
                                x = list(ncovs=1, xnames="treatment"),
                                X = list(treatment=-5),
                                nonprop = TRUE,
                                prior_sdnp = p_gamma(2,1),
                                tmin=0, tmax=10, nsim=10)
    p
  }, NA)
})

test_that("prior hazard SD over time", {
  expect_error({
    ## we allow the hazard to be extremely variable
    prior_haz_sd(mspline=sp,
                 coefs_mean = NULL,
                 prior_smooth = p_gamma(2,1),
                 prior_loghaz = p_normal(0,1),
                 nsim = 100, quantiles=c(0.25, 0.75))

    ## we are confident that the hazard will be close to constant
    prior_haz_sd(mspline=sp,
                 coefs_mean = NULL,
                 prior_smooth = p_gamma(1,40),
                 prior_loghaz = p_normal(0,1),
                 nsim = 100, quantiles=c(0.25, 0.75))

    ## proportional hazards model with one covariate, uncertainty in the cov effect
    ## Hazards are sometimes higher, then the SDs over time will be higher too.
    ## SD for hazard depends on the level of the hazard.
    prior_haz_sd(mspline=sp,
                 coefs_mean = NULL,
                 prior_smooth = p_gamma(1,40),
                 prior_loghaz = p_normal(0,1),
                 x = list(ncovs=1, xnames="treatment"), X = list(treatment=1),
                 nsim = 100, quantiles=c(0.25, 0.75))

    ## non proportional hazards model with one covariate. Even more uncertainty
    prior_haz_sd(mspline=sp,
                 coefs_mean = NULL,
                 prior_smooth = p_gamma(1,40),
                 prior_loghaz = p_normal(0,1),
                 x = list(ncovs=1, xnames="treatment"), X = list(treatment=1),
                 nonprop = TRUE,
                 prior_sdnp = p_gamma(2,1),
                 nsim = 100, quantiles=c(0.25, 0.75))

    plot_mspline(iknots = c(1,2,3), bknots = c(0, 5), df=5, coefs=c(1, rep(10, 5), 1))
    mspline_plotdata(iknots = c(1,2,3), bknots = c(0, 5), df=5, coefs=c(1, rep(10, 5), 1))

  }, NA)
})

test_that("prior hazard ratio SD over time", {

  ## Tiny nonproportionality effect: SDs for HR are close to zero.
  set.seed(1)
  psd <- prior_hr_sd(mspline=sp,
                     coefs_mean = NULL,
                     prior_smooth = p_gamma(1,400),
                     prior_loghaz = p_normal(0,1),
                     x = list(ncovs=1, xnames="treatment"),
                     X = list(treatment=1),
                     X0 = list(treatment=0),
                     nonprop = TRUE,
                     prior_sdnp = p_gamma(1, 100),
                     nsim = 100, quantiles=c(0.25, 0.75))
  expect_lt(psd$sd_hr[2], 0.5)

})

test_that("prior predictive functions returned with a fitted model",{
  mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
  X <- list(rxLev=1,"rxLev+5FU"=0)
  X0 <- list(rxLev=0,"rxLev+5FU"=0)
  expect_error({
    mod$prior_pred$haz(X = X, nsim=4)
    mod$prior_pred$haz_sd(X = X, quantiles = c(0.25, 0.75))
    mod$prior_pred$hr_sd(X = X, X0 = X0)
  }, NA)
})
