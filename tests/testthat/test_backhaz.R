
test_that("get_cum_backhaz",{
  backhaz <- data.frame(hazard = c(0.01, 0.02, 0.03), time = c(0, 5, 10))
  expect_equal(get_cum_backhaz(2, backhaz), 0.02)
  expect_equal(get_cum_backhaz(6, backhaz), 0.07)
  expect_equal(get_cum_backhaz(12, backhaz), 0.21)
  expect_equal(get_cum_backhaz(c(2,6,12), backhaz), c(0.02, 0.07, 0.21))
})


alpha <- -0.35
coefs <- c(0.0012, 0.029, 0.1, 0.15, 0.11, 0.11, 0.23, 0.0074, 0.052, 0.1)
knots <- c(0, 0.36, 0.62, 0.91, 1.2, 1.54, 2.02, 3)

test_that("distribution functions with backhaz",{
  backhaz <- data.frame(hazard = c(0.01, 0.02, 0.03), time = c(0, 5, 10))
  h1 <-  hsurvmspline(5, alpha, coefs, knots)
  h2 <- hsurvmspline(5, alpha, coefs, knots, backhaz=backhaz)
  expect_lt(h1, h2)
})
