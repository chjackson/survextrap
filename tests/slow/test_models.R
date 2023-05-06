postest_test <- function(x, newdata=NULL){
    niter <- 10
    t <- c(1, 5)
    summary(x)
    rmst(x, t=t, newdata=newdata, niter=niter)
    survival(x, newdata=newdata, t=t)
    hazard(x, newdata=newdata, t=t)
    invisible()
}

test_median <- function(mod, vname, value, tol=1e-01){
    expect_equal(summary(mod) %>% filter(variable==vname) %>% pull(median) %>% as.numeric(),
                 value, tol=tol)
}

nd <- data.frame(rx = c("Obs", "Lev+5FU"))

test_that("Basic spline model, no covariates",{
    ## VB fit. LOO warning
    suppressWarnings(modv <- survextrap(Surv(years, status) ~ 1,
                                        data=colons, fit_method="vb", loo=FALSE))
    expect_equal(coef(modv)["alpha"], coef(modm)["alpha"], tol=3e-01)
})
