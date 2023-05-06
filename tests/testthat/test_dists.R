alpha <- -0.35
coefs <- c(0.0012, 0.029, 0.1, 0.15, 0.23, 0.0074, 0.052, 0.1)
knots <- c(0.36, 0.62, 0.91, 1.2, 1.54, 2.02, 3)

test_that("psurvmspline",{
    q <- c(-1, 0, Inf, NA, NaN)
    expect_equal(psurvmspline(q, alpha, coefs, knots), c(0, 0, 1, NA, NaN))
    expect_equal(psurvmspline(q, alpha, coefs, knots, lower.tail = FALSE), c(1, 1, 0, NA, NaN))
    expect_equal(psurvmspline(q, alpha, coefs, knots, log.p = TRUE), c(-Inf, -Inf, 0, NA, NaN))
    q <- c(1.2, 1.3)
    pr <- psurvmspline(q, alpha, coefs, knots)
    qsurvmspline(pr, alpha, coefs, knots)
    expect_equal(qsurvmspline(pr, alpha, coefs, knots), q)

    q <- array(1, dim=c(2,3))
    pr <- psurvmspline(q, alpha, coefs, knots)
    expect_equivalent(attributes(q), attributes(pr))
})

test_that("qsurvmspline",{
    p <- c(0.1, 0.2)
    qu <- qsurvmspline(p, alpha, coefs, knots)
    expect_equal(psurvmspline(qu, alpha, coefs, knots), p)
})

test_that("dsurvmspline",{
    x <- c(-1, NA, NaN)
    expect_equal(dsurvmspline(x, alpha, coefs, knots), c(0, NA, NaN))
    x <- 1:3
    dens <- dsurvmspline(x, alpha, coefs, knots)
    surv <- psurvmspline(x, alpha, coefs, knots, lower.tail=FALSE)
    haz <- hsurvmspline(x, alpha, coefs, knots)
    expect_equal(haz, dens/surv)
})

test_that("hsurvmspline",{
    x <- c(-1, NA, NaN)
    expect_equal(hsurvmspline(x, alpha, coefs, knots), c(0, NA, NaN))
})

test_that("Hsurvmspline",{
    x <- c(-1, 0, Inf, NA, NaN)
    expect_equal(Hsurvmspline(x, alpha, coefs, knots), c(0, 0, NaN, NA, NaN))
    x <- 1:3
    cumhaz <- Hsurvmspline(x, alpha, coefs, knots)
    surv <- psurvmspline(x, alpha, coefs, knots, lower.tail = FALSE)
    expect_equal(surv, exp(-cumhaz))
})

test_that("rsurvmspline",{
    set.seed(1)
    ran <- rsurvmspline(10, alpha, coefs, knots)
    expect_type(ran, "double")
})

test_that("offsets in distribution functions",{
  q <- c(1.2)
  pr1 <- psurvmspline(q, alpha, coefs, knots)
  pr2 <- psurvmspline(q, alpha, coefs, knots, offsetH = 0.01)
  expect_lt(pr1, pr2)

  pr1 <- hsurvmspline(q, alpha, coefs, knots)
  pr2 <- hsurvmspline(q, alpha, coefs, knots, offseth = 0.01)
  expect_lt(pr1, pr2)

  pr1 <- Hsurvmspline(q, alpha, coefs, knots)
  pr2 <- Hsurvmspline(q, alpha, coefs, knots, offsetH = 0.01)
  expect_lt(pr1, pr2)

  pr1 <- dsurvmspline(q, alpha, coefs, knots)
  pr2 <- dsurvmspline(q, alpha, coefs, knots, offseth=0.01, offsetH = 0.02)
  expect_true(pr1 != pr2)

  q <- c(1.2, 1.3)
  pr1 <- psurvmspline(q, alpha, coefs, knots, offsetH = c(NA, 0.01))
  pr2 <- psurvmspline(1.3, alpha, coefs, knots, offsetH = c(0.01))
  expect_equal(pr1[2], pr2)

  pr1 <- dsurvmspline(1.2, alpha, coefs, knots, offseth=0.01, offsetH=0.02)
  pr2 <- dsurvmspline(1.2, alpha, coefs, knots, offseth=0.01, offsetH = c(NA,0.02))
  expect_equal(pr1, pr2[2])
})

test_that("dists with constant hazard",{
    knots <- c(1,2,3)
    cf <- mspline_constant_coefs(list(knots=knots))
    haz <- hsurvmspline(c(1), 0, cf, knots)
    cumhaz <- Hsurvmspline(3, 0, cf, knots)
    expect_equal(haz * 3, cumhaz)
    prob <- psurvmspline(3, 0, cf, knots, lower.tail = FALSE)
    expect_equal(prob, exp(-cumhaz))

    alpha <- c(0, 1)
    cf2 <- rbind(cf, cf)
    haz <- hsurvmspline(c(1), alpha, cf2, knots)
    cumhaz <- Hsurvmspline(3, alpha, cf2, knots)
    expect_equal(haz * 3, cumhaz)
    prob <- psurvmspline(3, alpha, cf2, knots, lower.tail = FALSE)
    expect_equal(prob, exp(-cumhaz))
})
