#' Simulate a dataset from the prior predictive distribution of survival times 
#' in an M-spline survival model.
#'
#' Simulate a dataset from the prior predictive distribution of
#' survival times in an M-spline survival model.  Additive hazards
#' models not currently supported.
#'
#' @param n Sample size of the simulated dataset.  Each observation in
#'   the dataset is generated from a model with the same parameters.
#'   These parameters are generated from a single simulation from the
#'   prior distribution.
#'
#' @param fix_prior If \code{TRUE}, then one value of the parameter
#'   vector is drawn from the prior, followed by \code{n}
#'   individual-level times given this common prior value.  If
#'   \code{FALSE}, then to produce each sampled individual time, a
#'   different sample from the prior is used.
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
                       fix_prior=FALSE,
                       mspline,
                       censtime = Inf,
                       coefs_mean = NULL,
                       prior_hscale = p_normal(0,20),
                       prior_hsd = p_gamma(2,1),
                       X = NULL,
                       prior_loghr = NULL,
                       prior_hrsd = NULL,
                       prior_cure = NULL)
{
  mspline <- mspline_default(mspline)
  nprior <- if (fix_prior) 1 else n
  sam <- prior_sample(mspline = mspline,
                      coefs_mean = coefs_mean,
                      prior_hscale = prior_hscale,
                      prior_hsd = prior_hsd,
                      prior_loghr = prior_loghr,
                      prior_hrsd = prior_hrsd,
                      prior_cure = prior_cure,
                      X = X, X0=NULL, nsim=nprior)
  sim <- rsurvmspline(n, alpha=sam$alpha, coefs=sam$coefs,
                      knots=mspline$knots, degree=mspline$degree,
                      pcure=sam$pcure,
                      bsmooth=mspline$bsmooth) # TODO backhaz 
  res <- data.frame(time = pmin(sim, censtime),
                    event = as.numeric(sim<=censtime))
  attr(res, "prior") <- sam
  res
}	
