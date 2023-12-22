
mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")

test_that("Default newdata for factor covariates",{
  mn <- mean(mod, niter=10)
  expect_equal(nrow(mn), length(levels(colons$rx)))
})

test_that("Errors in specifying newdata",{
  expect_error(mean(mod, niter=10, newdata=1),
               "`newdata` should be a data frame")
  ndwrong <- data.frame(blah = 1)
  expect_error(mean(mod, newdata=ndwrong),
               "Values of covariate `rx` not included in `newdata`")
  ndwrong <- data.frame(blah = 1, rx="Wrong level")
  expect_error(mean(mod, newdata=ndwrong),
               "factor rx has new level")
})

test_that("Different summary functions",{
  h1 <- hazard(mod, niter=10, t=1)
  h2 <- hazard(mod, niter=10, t=1, summ_fns = list(mean=mean, median=median))
  expect_equal(h1$median, h2$median)
  h2 <- hazard(mod, niter=10, t=1,
               summ_fns = list(mean=mean, ~quantile(.x, c(0.2, 0.8))))
  expect_lt(h2$`20%`[1], h2$`80%`[1])
  h2 <- hazard(mod, niter=10, t=1,
               summ_fns = list(mean=mean, ~quantile(.x, c(0.025, 0.975))))
  expect_equal(h2$`2.5%`, h1$lower)
  s1 <- survival(mod, niter=10, t=1)
  s2 <- survival(mod, niter=10, t=1,
               summ_fns = list(mean=mean, ~quantile(.x, c(0.025, 0.975))))
  expect_equal(s2$`2.5%`, s1$lower)
  r1 <- rmst(mod, t=5, niter=10)
  r2 <- rmst(mod, t=5, niter=10,
             summ_fns = list(mean=mean, ~quantile(.x, c(0.025, 0.975))))
  expect_equal(r1$lower, r2$`2.5%`)
})


nd <- data.frame(rx = c("Lev+5FU","Lev"))

test_that("Hazard ratio", {
  hr <- hazard_ratio(mod, newdata=nd, t=1:2, niter=20)
  expect_equal(hr$median[1], hr$median[2])
  plot_hazard_ratio(mod, newdata=nd, niter=20, t=1:5)
  hr2 <- hazard_ratio(mod, newdata=nd, t=1:2, niter=20,
                      summ_fns = list(mean=mean))
  hr3 <- hazard_ratio(mod, newdata=nd, t=1:2, niter=20,
                      summ_fns = list(mean=mean, median=median))
  expect_equal(hr2$mean, hr3$mean)
  expect_equal(hr$median, hr3$median)
  hr4 <- hazard_ratio(mod, newdata=nd, t=1:2, niter=20,
                      summ_fns = list(~quantile(.x, c(0.2, 0.8))))
  expect_lt(hr4$`20%`[1], hr4$`80%`[1])
})

test_that("hrtime", {
  hrt <- hrtime(mod, newdata=nd, niter=30)
  expect_equal(hrt[1,"median"], hrt[2,"median"]) # because proportional hazards
  hrt2 <- hrtime(mod, newdata=nd, niter=30,
                summ_fns = list(mean=mean, median=median))
  expect_equal(hrt2$median, hrt$median)
  hrt3 <- hrtime(mod, newdata=nd[1,,drop=FALSE], niter=30,
                 summ_fns = list(mean=mean))
  expect_equal(hrt3$mean, hrt2$mean[1])
})
