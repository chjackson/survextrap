#' Simulate a dataset from the prior predictive distribution of survival times 
#' in an M-spline survival model.
#'
#' @param n Sample size of the simulated dataset.  Each observation in the dataset is
#' generated from a model with the same parameters.  These parameters are generated from
#' a single simulation from the prior distribution.
#'
#' @param censtime Right-censoring time to impose on the simulated event times.
#'
#' @param alpha_mean Prior mean of the normal distribution for the baseline log hazard.
#'
#' @param alpha_sd Prior standard deviation of the normal distribution for the baseline log hazard.
#'
#' @param lcoef_mean Prior mean vector of the the logistic distribution for the spline coefficients on the multinomial logit scale.
#' 
#' @param lcoef_sd Prior standard deviation of the logistic distribution for the spline coefficients on the multinomial logit scale.
#'
#' @param iknots Internal knots of the M-spline.
#'
#' @param bknots Boundary knots of the M-spline.
#'
#' @param degree M-spline polynomial degree.
#'
#' @return A data frame with columns \code{time} (simulated time) and \code{event} (indicator for whether
#' the time is an event time, as opposed to a right-censoring time).  The prior parameters are returned
#' in the \code{prior} attribute as a list with components \code{alpha} (baseline log hazard) and
#' \code{coefs} (spline coefficients).
#'
#' @export
prior_pred <- function(n, censtime=Inf, 
                       alpha_mean=0, alpha_sd=20, lcoef_mean = NULL, lcoef_sd=1,
                       iknots, bknots, degree=3){
  nsim <- 1
  if (is.null(lcoef_mean))
    lcoef_mean <- mspline_uniform_weights(iknots, bknots, logit=TRUE)
  np <- length(lcoef_mean) + 1
  beta <- matrix(nrow=nsim, ncol=np)
  beta[,1] <- 0
  alpha <- exp(rnorm(nsim, alpha_mean, alpha_sd))
  for (j in 2:np){
    beta[,j] <- rlogis(nsim, lcoef_mean[j-1], lcoef_sd)
  }
  coefs <- exp(beta) / rowSums(exp(beta))
  knots <- c(iknots, bknots)
  sim <- rsurvmspline(n, alpha, coefs, knots, degree=degree)
  res <- data.frame(time = pmin(sim, censtime), event = as.numeric(sim<=censtime))
  attr(res, "prior") <- list(alpha=alpha, coefs=coefs)
  res
}	

