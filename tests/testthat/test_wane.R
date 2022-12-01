alpha <- -0.35
coefs <- c(0.0012, 0.029, 0.1, 0.15, 0.11, 0.11, 0.23, 0.0074, 0.052, 0.1)
knots <- c(0, 0.36, 0.62, 0.91, 1.2, 1.54, 2.02, 3)

psurvmspline(q=10, alpha=alpha, coefs=coefs, knots=knots)
psurvmspline_wane(q=10, alpha1=alpha, alpha0=-0.2, coefs=coefs, knots=knots, wane_period=c(11, 15))
psurvmspline_wane(q=10, alpha1=alpha, alpha0=alpha, coefs=coefs, knots=knots, wane_period=c(9, 11))
psurvmspline_wane(q=10, alpha1=alpha, alpha0=alpha*1.01, coefs=coefs, knots=knots, wane_period=c(9, 11))

# vectorised
psurvmspline_wane(q=c(0, 4, 8, 10, 12), alpha1=alpha, alpha0=alpha*1.01, coefs=coefs, knots=knots,
                  wane_period=c(9, 11))

Hsurvmspline_wane(x=c(0, 8, 10, 12), alpha1=alpha, alpha0=alpha*1.01, coefs=coefs, knots=knots,
                  wane_period=c(9, 11))

hsurvmspline_wane(x=c(0, seq(8, 12, by=0.1)), alpha1=alpha, alpha0=-5, coefs=coefs, knots=knots,
                  wane_period=c(9, 11))

# vectorisation over parameters

expect_error(Hsurvmspline_wane(x=c(0, 8, 10, 12), alpha1=rep(alpha,2), alpha0=alpha*1.01,
                               coefs=coefs, knots=knots, wane_period=c(9, 11)),
             "lengths of `alpha0`")

h1 <- Hsurvmspline_wane(x=c(8, 10), alpha1=-1, alpha0=-2,
                        coefs=coefs, knots=knots, wane_period=c(9, 11))
h2 <- Hsurvmspline_wane(x=c(8, 10), alpha1=-3, alpha0=-4,
                        coefs=coefs, knots=knots, wane_period=c(9, 11))
h12 <- Hsurvmspline_wane(x=c(8, 10), alpha1=c(-1,-3), alpha0=c(-2,-4),
                         coefs=coefs, knots=knots, wane_period=c(9, 11))
expect_equal(h12, c(h1[1], h2[2]))

