#' Simulate a dataset from the prior predictive distribution of survival times 
#' in an M-spline survival model.
#'
#' Simulate a dataset from the prior predictive distribution of survival times 
#' in an M-spline survival model.  Mixture cure and additive hazards models
#' not currently supported.
#'
#' @param n Sample size of the simulated dataset.  Each observation in
#'   the dataset is generated from a model with the same parameters.
#'   These parameters are generated from a single simulation from the
#'   prior distribution.
#'
#' @param censtime Right-censoring time to impose on the simulated
#'   event times.
#'
#' @inheritParams prior_sample
#'
#' @return A data frame with columns \code{time} (simulated time) and
#'   \code{event} (indicator for whether the time is an event time, as
#'   opposed to a right-censoring time).  The prior parameters are
#'   returned in the \code{prior} attribute as a list with components
#'   \code{alpha} (baseline log hazard) and \code{coefs} (spline
#'   coefficients).
#'
#' @seealso \code{\link{prior_sample}}
#'
#' @export
prior_pred <- function(n,
                       mspline,
                       censtime = Inf,
                       coefs_mean = NULL,
                       prior_hsd = p_gamma(2,1),
                       prior_hscale,
                       X = NULL,
                       prior_loghr = NULL,
                       prior_hrsd = NULL,
                       prior_cure = NULL)
{
  sam <- prior_sample(mspline = mspline,
                      coefs_mean = coefs_mean,
                      prior_hsd = prior_hsd,
                      prior_hscale = prior_hscale,
                      prior_hrsd = prior_hrsd,
                      X = X, X0=NULL, nsim=n)
  sim <- rsurvmspline(n, alpha=sam$alpha, coefs=sam$coefs,
                      knots=mspline$knots, degree=mspline$degree,
                      bsmooth=mspline$bsmooth) # TODO pcure.  backhaz too? 
  res <- data.frame(time = pmin(sim, censtime),
                    event = as.numeric(sim<=censtime))
  attr(res, "prior") <- sam
  res
}	
