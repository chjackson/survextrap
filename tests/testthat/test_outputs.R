
mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")

test_that("Default newdata for factor covariates",{
  mn <- mean(mod, niter=10)
  expect_equal(nrow(mn), length(levels(colons$rx)))
})

test_that("Errors in specifying newdata",{
  expect_error(mean(mod, niter=10, newdata=list(rx=c("Obs","Lev"))),
               "`newdata` should be a data frame")
  expect_error(mean(mod, niter=10, newdata=1),
               "`newdata` should be a data frame")
  ndwrong <- data.frame(blah = 1)
  expect_error(mean(mod, newdata=ndwrong),
               "Values of covariate `rx` not included in `newdata`")
  ndwrong <- data.frame(blah = 1, rx="Wrong level")
  expect_error(mean(mod, newdata=ndwrong),
               "factor rx has new level")
})
