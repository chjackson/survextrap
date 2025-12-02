
test_that("spline demo plots", {
  expect_error(mspline_plotdata(knots = c(1,2,3,5), df=5, bsmooth=FALSE,
                                coefs=c(1, rep(10, 3), 1)),
               "dimension of `coefs` is 5, should be 7")
})

sp <- mspline_init(knots=c(1,2,3,5), degree=3)

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
                 formula = ~ treatment,
                 newdata = list(treatment=1),
                 nsim = 100)

    ## 2 covs
    sam1 <- prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_loghr = list(treatment=p_normal(0,1),
                                    age=p_normal(1,2)),
                 newdata = list(treatment=1, age=1),
                 formula = ~treatment + age,
                 nsim = 100)

    ## newdata with more than one value is OK
    sam <-  prior_sample(mspline=sp,
                         coefs_mean = NULL,
                         prior_hsd = p_gamma(2,1),
                         prior_hscale = p_normal(0,1),
                         prior_loghr = list(treatment=p_normal(0,1),
                                            age=p_normal(1,2)),
                         newdata = list(treatment=c(0,1), age=1:2),
                         formula = ~treatment + age,
                         nsim = 100)

    ## .. with an interaction
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_loghr = list(treatment=p_normal(0,1),
                                    age=p_normal(1,2),
                                    "treatment:age"=p_normal(0,0.1)),
                 newdata = list(treatment=1, age=1),
                 formula = ~treatment*age,
                 nsim = 100)

    expect_error(prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_loghr = list(age=p_normal(0,1)),
                 newdata = list(treatment=1),
                 formula = ~treatment,
                 nsim = 100), "names of prior_loghr do not match")

    expect_error(prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_loghr = list(treatment=p_normal(1, 2), age=p_normal(0,1)),
                 newdata = list(treatment=1),
                 formula = ~treatment,
                 nsim = 100), "loghr has 2 components, but there is 1 coefficient")

    expect_error(prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_loghr = list(treatment=p_normal(1, 2)),
                 newdata = list(treatment=1),
                 nsim = 100), "regression formula should be supplied in `formula`")

    expect_error(prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_loghr = list(treatment=p_normal(1, 2)),
                 formula = ~treatment,
                 nsim = 100), "covariate values should be supplied in `newdata`")

    expect_error(prior_sample(mspline=sp,
                              coefs_mean = NULL,
                              prior_hsd = p_gamma(2,1),
                              prior_hscale = p_normal(0,1),
                              prior_loghr = list(treatment=p_normal(1, 2)),
                              nsim = 100), "regression formula should be supplied in `formula`")

    ## nonprop haz model, 1 cov
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_hrsd = p_gamma(2,1),
                 newdata = list(treatment=1),
                 formula = ~treatment,
                 nonprop = ~treatment,
                 nsim = 100)

    ## nonprop haz model, 2 covs
    prior_sample(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_hrsd = p_gamma(2,1),
                 newdata = list(treatment=1),
                 formula = ~treatment,
                 nonprop = ~treatment,
                 nsim = 100)

    ## cure, no covs
    prior_sample(mspline = sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_cure = p_beta(2,10),
                 nsim = 100)

    ## errors with cure and covs
    expect_error(
    prior_sample(mspline = sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 prior_cure = p_beta(2,10),
                 prior_logor_cure = p_normal(0,2),
                 nsim = 100), "a regression formula should be supplied in `cure`")

    expect_error(
      prior_sample(mspline = sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 cure = ~treatment,
                 prior_cure = p_beta(2,10),
                 prior_logor_cure = p_normal(0,2),
                 nsim = 100), "covariate values should be supplied"
    )

    ## cure with covs
    prior_sample(mspline = sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(2,1),
                 prior_hscale = p_normal(0,1),
                 cure = ~treatment,
                 newdata = list(treatment = 1),
                 prior_cure = p_beta(2,10),
                 prior_logor_cure = p_normal(0,2),
                 nsim = 100)

  },NA)
})


test_that("prior hazard function samples", {
  expect_error({
    pdf <- prior_sample_hazard(knots = c(1,2,3,5), degree=3,
                               prior_hscale = p_normal(0,1),
                               prior_hsd = p_gamma(10, 10),
                               tmin=0, tmax=10, nsim=10)
    ggplot(pdf, aes(x=time, y=haz, group=rep)) + geom_line() + ylim(0,2)

    pdf <- prior_sample_hazard(knots = c(1,2,3,5), degree=3,
                               prior_hscale = p_normal(0,1),
                               prior_hsd = p_gamma(10, 10),
                               smooth_model = "random_walk",
                               tmin=0, tmax=10, nsim=10)
    ggplot(pdf, aes(x=time, y=haz, group=rep)) + geom_line() + ylim(0,2)

    ## prop haz model
    pdf <- prior_sample_hazard(knots = c(1,2,3,5), degree=3,
                               prior_hscale = p_normal(0,1),
                               prior_hsd = p_gamma(10, 10),
                               formula = ~ treatment,
                               newdata = list(treatment=-5),
                               tmin=0, tmax=10, nsim=10)
    ggplot(pdf, aes(x=time, y=haz, group=rep)) + geom_line() + ylim(0,2)

    ## nonprop haz model
    p <- plot_prior_hazard(knots = c(1,2,3,5), degree=3,
                           prior_hscale = p_normal(0,1),
                           prior_hsd = p_gamma(10, 10),
                           formula = ~treatment,
                           newdata = list(treatment=-5),
                           nonprop = ~treatment,
                           prior_hrsd = p_gamma(2,1),
                           tmin=0, tmax=10, nsim=10)
    p

    ## multiple covariates with different priors - must be named
    psam <- prior_sample(mspline = mspline_init(knots = c(1,2,3,5), degree=3),
                         prior_hscale = p_normal(0,1),
                         prior_hsd = p_gamma(10, 10),
                         prior_loghr = list(age=p_normal(2,3),
                                            treatment=p_normal(0,0.01)),
                         formula = ~treatment + age,
                         newdata = list(treatment=-5, age=50),
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
                 formula = ~ treatment,
                 newdata = list(treatment=1),
                 nsim = 100, quantiles=c(0.25, 0.75))

    ## non proportional hazards model with one covariate. Even more uncertainty
    prior_haz_sd(mspline=sp,
                 coefs_mean = NULL,
                 prior_hsd = p_gamma(1,40),
                 prior_hscale = p_normal(0,1),
                 formula = ~ treatment,
                 newdata = list(treatment=1),
                 nonprop = ~treatment,
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
                     newdata = list(treatment=1),
                     newdata0 = list(treatment=0),
                     formula = ~treatment,
                     nonprop = ~treatment,
                     prior_hrsd = p_gamma(1, 100),
                     nsim = 100, quantiles=c(0.25, 0.75))
  expect_lt(psd$sd_hr[2], 0.5)

})

test_that("prior sampling functions returned with a fitted model",{
  mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
  levs <- c("Obs","Lev","Lev+5FU")
  newdata <- data.frame(rx=factor("Lev", levels = levs))
  newdata0 <- data.frame(rx=factor("Lev+5FU", levels = levs))
  expect_error({
    mod$prior_sample$haz(newdata=newdata, nsim=4)
    mod$prior_sample$haz_sd(newdata=newdata, quantiles = c(0.25, 0.75))
    mod$prior_sample$hr_sd(newdata=newdata, newdata0=newdata0)
  }, NA)
  mod <- survextrap(Surv(years, status) ~ rx, data=colons, cure=TRUE, fit_method="opt")
  expect_error({
    mod$prior_sample$haz_sd(newdata=newdata, nsim=4)
  }, NA)
  mod <- survextrap(Surv(years, status) ~ rx, data=colons, cure=~rx, fit_method="opt")
  expect_error({
    mod$prior_sample$haz_sd(newdata=newdata, nsim=4)
  }, NA)
})

test_that("prior_haz_const interprets hazard scale prior",{
  p1 <- prior_haz_const(sp)
  p2 <- prior_haz_const(sp, prior_hscale = p_normal(0,1))
  expect_lt(p1$haz[1], p2$haz[1])
})
