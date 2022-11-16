
sim_rank <- function(i, fit_method="opt"){
  library(survextrap) # have to do this to parallelise on Windows, inefficient
  L <- 100
  iknots <- 1:3; bknots=c(0,5); knots <- c(iknots, bknots)
  simdf <- prior_pred(n=200, censtime=5, 
                      alpha_mean=0, alpha_sd=1, lcoef_mean = NULL, lcoef_sd=1, 
                      iknots=iknots, bknots=bknots, degree=3)
  iter <- if (fit_method=="mcmc") 1000 else 2000
  open_progress <- if(fit_method=="mcmc") FALSE else NULL
  tryres <- try(mod <- survextrap(Surv(time, event) ~ 1, data=simdf, fit_method=fit_method,
                                  iter = iter, 
                                  prior_loghaz = p_normal(0, 20), 
                                  basehaz_ops=list(iknots=iknots,bknots=bknots,df=7),
                                  loo = FALSE,
                                  refresh = 0))
  if (!inherits(tryres, "try-error")){
    summ <- attr(rmst(mod, t=10, niter=L),"sample") # TODO ensure at least L samples
    prior <- attr(simdf, "prior")
    summ_prior <- rmst_survmspline(t=10, prior$alpha, prior$coefs, knots)
    res <- findInterval(summ_prior, sort(summ))
  } else res <- NA
  res
}
