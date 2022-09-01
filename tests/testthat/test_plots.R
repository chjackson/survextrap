mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
nd <- data.frame(rx = c("Lev+5FU","Lev"))

test_that("hazard ratio", {
  hr <- hazard_ratio(mod, newdata=nd, times=1:2, niter=20)
  expect_equal(hr$median[1], hr$median[2])
  plot_hazard_ratio(mod, newdata=nd, niter=20, times=1:5)
})

test_that("survival and hazard plots", {
  expect_error({
    p <- plot_survival(mod, niter=20)
    p <- plot_survival(mod, niter=20, km=FALSE)
    p <- plot_hazard(mod, niter=20)
    p <- plot(mod, niter=20)
  }, NA)
})

test_that("deconstruct fitted spline",{
  expect_error({
    hazdf <-  deconstruct_mspline(mod)
    bknots <- unname(attr(hazdf,"bknots"))
    iknots <- unname(attr(hazdf,"iknots"))
    p <- ggplot(hazdf, aes(x=time, y=value, group=term)) +
      geom_line(alpha=0.5) +
      scale_x_continuous(breaks=c(iknots, bknots)) +
      theme_minimal() +
      theme(panel.grid.minor = element_blank()) +
      geom_vline(xintercept = bknots, col="gray50") +
      xlab("Time") + ylab("Hazard")
  }, NA)
})

test_that("spline demo plots", {
  expect_error(mspline_plotdata(iknots = c(1,2,3), bknots = c(0, 5), df=5, coefs=c(1, rep(10, 3), 1)),
               "length of `coefs` is 5, should be 7")

  expect_error({
    plot_mspline(iknots = c(1,2,3), bknots = c(0, 5), df=5, coefs=c(1, rep(10, 5), 1))
    mspline_plotdata(iknots = c(1,2,3), bknots = c(0, 5), df=5, coefs=c(1, rep(10, 5), 1))

    mspline_priorpred_df(iknots=1:3, bknots=c(0,5), df=5, degree=2,
                         prior_mean=c(1, rep(10,3), 1), prior_sd=1, scale=1, scale_sd=1,
                         tmin=0, tmax=10, nsim=10)
    p <- plot_mspline_priorpred(iknots=1:3, bknots=c(0,5), df=5, degree=2,
                                prior_mean=c(1, rep(10,3), 1), prior_sd=1, scale=1, scale_sd=1,
                                tmin=0, tmax=10, nsim=10)
  }, NA)
})
