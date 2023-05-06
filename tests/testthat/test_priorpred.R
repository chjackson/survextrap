sp <- list(knots=c(1,2,3,5), degree=3)

test_that("prior prediction", {
  set.seed(1)
  pred <- prior_pred(10, mspline = sp,
                     coefs_mean = NULL,
                     prior_hsd = p_gamma(2,1),
                     prior_hscale = p_normal(0,1),
                     censtime = 20)
  expect_is(pred, "data.frame")
})
