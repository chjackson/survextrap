
test_that("spline demo plots", {
  expect_error(mspline_plotdata(knots = c(1,2,3,5), df=5, bsmooth=FALSE,
                                coefs=c(1, rep(10, 3), 1)),
               "dimension of `coefs` is 5, should be 7")
})

sp <- list(knots=c(1,2,3,5), degree=3)

test_that("prior hazard parameter samples", {
  expect_error({ # ?? any sensible expectations here?
    ## no covariates
    prior_sample(mspline = sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 nsim = 100)

    ## prop haz model, 1 cov
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 X = list(treatment=1),
                 nsim = 100)
    ## prop haz model, 2 covs
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 X = list(treatment=1, age=60),
                 nsim = 100)

    ## nonprop haz model, 1 cov
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_hrsd = p_gamma(2,1),
                 X = list(treatment=1),
                 nsim = 100)
    ## nonprop haz model, 2 covs
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_hrsd = p_gamma(2,1),
                 X = list(treatment=1),
                 nsim = 100)

    prior_sample(mspline = sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_cure = p_beta(2,10),
                 nsim = 100)

  },NA)
})


test_that("prior hazard function samples", {
  expect_error({
    pdf <- mspline_priorpred(knots = c(1,2,3,5), degree=3,
                                prior_hscale = p_normal(0,1),
                                prior_hsd = p_gamma(10, 10),
                                tmin=0, tmax=10, nsim=10)
    ggplot(pdf, aes(x=time, y=haz, group=rep)) + geom_line() + ylim(0,2)

    ## prop haz model
    pdf <- mspline_priorpred(knots = c(1,2,3,5), degree=3,
                                prior_hscale = p_normal(0,1),
                                prior_hsd = p_gamma(10, 10),
                                X = list(treatment=-5),
                                tmin=0, tmax=10, nsim=10)
    ggplot(pdf, aes(x=time, y=haz, group=rep)) + geom_line() + ylim(0,2)

    ## nonprop haz model
    p <- plot_mspline_priorpred(knots = c(1,2,3,5), degree=3,
                                prior_hscale = p_normal(0,1),
                                prior_hsd = p_gamma(10, 10),
                                X = list(treatment=-5),
                                prior_hrsd = p_gamma(2,1),
                                tmin=0, tmax=10, nsim=10)
    p

    ## multiple covariates with different priors - must be named
    psam <- prior_sample(mspline = list(knots = c(1,2,3,5), degree=3),
                         prior_hscale = p_normal(0,1),
                         prior_hsd = p_gamma(10, 10),
                         prior_loghr = list(age=p_normal(2,3),
                                            treatment=p_normal(0,0.01)),
                         X = list(treatment=-5, age=50),
                         nsim=10)

  }, NA)

})

test_that("prior hazard SD over time", {
  expect_error({
    ## we allow the hazard to be extremely variable
    prior_haz_sd(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 nsim = 100, quantiles=c(0.25, 0.75))

    ## we are confident that the hazard will be close to constant
    prior_haz_sd(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(1,40),
                 prior_hscale = p_normal(0,1),
                 nsim = 100, quantiles=c(0.25, 0.75))

    ## proportional hazards model with one covariate, uncertainty in the cov effect
    ## Hazards are sometimes higher, then the SDs over time will be higher too.
    ## SD for hazard depends on the level of the hazard.
    prior_haz_sd(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(1,40),
                 prior_hscale = p_normal(0,1),
                 X = list(treatment=1),
                 nsim = 100, quantiles=c(0.25, 0.75))

    ## non proportional hazards model with one covariate. Even more uncertainty
    prior_haz_sd(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(1,40),
                 prior_hscale = p_normal(0,1),
                 X = list(treatment=1),
                 prior_hrsd = p_gamma(2,1),
                 nsim = 100, quantiles=c(0.25, 0.75))

    plot_mspline(knots = c(1,2,3,5), coefs=c(1, rep(10, 5), 1), bsmooth = FALSE)
    mspline_plotdata(knots = c(1,2,3,5), coefs=c(1, rep(10, 5), 1), bsmooth = FALSE)

  }, NA)
})

test_that("prior hazard ratio SD over time", {

  ## Tiny nonproportionality effect: SDs for HR are close to zero.
  set.seed(1)
  psd <- prior_hr_sd(mspline=sp,
                     coefs_mean = NULL,
                     prior_hsd = p_gamma(1,400),
                     prior_hscale = p_normal(0,1),
                     X = list(treatment=1),
                     X0 = list(treatment=0),
                     prior_hrsd = p_gamma(1, 100),
                     nsim = 100, quantiles=c(0.25, 0.75))
  expect_lt(psd$sd_hr[2], 0.5)

})

test_that("prior sampling functions returned with a fitted model",{
  mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
  X <- list(rxLev=1,"rxLev+5FU"=0)
  X0 <- list(rxLev=0,"rxLev+5FU"=0)
  expect_error({
    mod$prior_sample$haz(X = X, nsim=4)
    mod$prior_sample$haz_sd(X = X, quantiles = c(0.25, 0.75))
    mod$prior_sample$hr_sd(X = X, X0 = X0)
  }, NA)
})
