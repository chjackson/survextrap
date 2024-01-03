#' Constructor for a standardising population used for survextrap
#' outputs
#'
#' Standardised outputs are outputs from models with covariates, that
#' are defined by marginalising (averaging) over covariate values in a
#' given population, rather than being conditional on a given
#' covariate value.
#'
#' @details These are produced by generating a Monte Carlo sample from
#'   the joint distribution of parameters \eqn{\theta} and covariate
#'   values \eqn{X}, \eqn{p(X,\theta) = p(\theta|X)p(X)}, where
#'   \eqn{p(X)} is defined by the empirical distribution of covariates
#'   in the standard population.
#'
#' Hence applying a vectorised output function \eqn{g()} (such as the
#' RMST or survival probability) to this sample produces a sample from
#' the posterior of \eqn{\int g(\theta|X) dX}: the average RMST (say)
#' for a heterogeneous population.
#'
#' @aliases standardize_to
#'
#' @details See the Examples vignette for some examples and notes on
#'   computation.
#' 
#' @param newdata Data frame describing a population.
#'
#' @param random By default this is \code{FALSE}, indicating that
#'   standardised samples should be obtained by concatenating the
#'   posterior samples for each covariate value in the standard
#'   population.  The sample from the standardised posterior of
#'   parameters then has size \code{niter} times the number of rows in
#'   \code{newdata}, where \code{niter} is the number of MCMC
#'   iterations used in the original \code{survextrap} fit.  Computing
#'   the resulting output function (e.g. RMST which uses numerical
#'   integration) can then be computationally intensive if this sample
#'   size is large.
#'
#'   A quicker alternative is to sample a random row of the standard
#'   population for each MCMC iteration.  The standardised sample from
#'   the posterior then has size \code{niter}.  This is specified by
#'   using \code{random=TRUE}.  If this is used, then the result
#'   depends on the random number seed, and it should be checked that
#'   the results are stable to within the required number of
#'   significant figures.  If not, run \code{survextrap} with more
#'   MCMC iterations or increase \code{nstd} here.
#'
#' @param nstd Number of draws from the population distribution used
#'   per MCMC sample from the parameters when \code{random=TRUE}.
#'   With the default of 1, the value of the covariate vector \eqn{X}
#'   is essentially treated as if it were an additional parameter in
#'   the Bayesian model, drawn by Monte Carlo independently of the
#'   remaining parameters.
#'
#' @return A copy of \code{newdata}, but with attributes added to
#'   indicate that this should be used as a standard population.  When
#'   this \code{newdata} is passed to \code{survextrap}'s output
#'   functions, the outputs will then be presented as an average over
#'   the empirical distribution of covariate values described by
#'   \code{newdata}, rather than as one output per row of
#'   \code{newdata} (distinct covariate values).
#'
#' @examples
#' rxph_mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
#' ref_pop <- data.frame(rx = c("Obs","Lev+5FU"))
#' 
#' # covariate-specific outputs
#' survival(rxph_mod, t = c(5,10), newdata = ref_pop)
#' 
#' # standardised outputs
#' survival(rxph_mod, t = c(5,10), newdata = standardise_to(ref_pop))
#'
#' @export
standardise_to <- function(newdata, nstd=1, random=FALSE){
  attr(newdata, "std") <- TRUE
  if (!is.numeric(nstd)) stop("`nstd` must be numeric")
  if (nstd < 1) stop("`nstd` must be 1 or more")
  attr(newdata, "nstd") <- nstd
  attr(newdata, "random") <- random
  newdata
}

#' @rdname standardise_to
#' @export
standardize_to <- standardise_to


#' Marginalise output of get_pars over a standard population
#'
#' Without marginalisation, \code{get_pars} returns, for each
#' parameter \eqn{\theta} and covariate value \eqn{X}, samples from
#' the posterior of \eqn{\theta|X}.  For each parameter, \code{nvals}
#' sample sets of size \code{niter} are produced, one for each set of
#' covariate values.  \code{niter} is the number of MCMC samples
#' requested, and \code{nvals} is the number of rows of
#' \code{newdata}.
#'
#' With marginalisation, \code{get_pars} returns a single set of
#' parameter values, sampled from the joint distribution of
#' \eqn{\theta|X} and \eqn{X}.
#'
#' Hence applying a vectorised function \eqn{g()} to the marginalised
#' sample produces a sample from the posterior of \eqn{\int
#' g(\theta|X) dX}.
#'
#' The distribution of \eqn{X} is defined by the empirical
#' distribution of values in the \code{newdata} supplied to
#' \code{get_pars}.
#'
#' @inheritParams standardise_to
#'
#' @param pars output of \code{get_pars}
#'
#' @return A list with components:
#'
#' \code{alpha} Log hazard scale.  Matrix with \code{niter} rows and
#' one column.
#'
#' \code{coefs} Spline basis coefficients.  Array with dimensions
#' \code{niter} x \code{1} x \code{nvars}, where \code{nvars} is the
#' basis dimension.
#'
#' \code{pcure} Cure probability in cure models. Matrix with
#' \code{niter} rows and one column.
#'
#' @noRd
marginalise_pars <- function(pars, nstd, random=FALSE){
  niter <- pars$niter
  nvals <- ncol(pars$alpha)
  nvars <- dim(pars$coefs)[3]

  if (random){
    alpha <- numeric(niter*nstd)
    coefs <- matrix(nrow=niter*nstd, ncol=nvars)
    for (i in 1:nstd){
      popind <- sample(1:nvals, niter, replace=TRUE)
      matind <- niter*(i-1) + 1:niter
      alpha[matind] <- pars$alpha[cbind(1:niter, popind)]
      inds3 <- cbind(rep(matind, nvars), popind, rep(1:nvars, each=niter))
      coefs[matind,] <- pars$coefs[inds3]
    }
    alpha_std <- matrix(alpha, ncol=1)
    coefs_std <- array(coefs, dim = c(niter, 1, nvars))
    if (!is.null(pars$pcure)){
      pcure <- numeric(niter*nstd)
      for (i in 1:nstd){
        pcure[matind] <- pars$pcure[cbind(1:niter, popind)]
      }
      pcure_std <- matrix(pcure, ncol=1)
    } else pcure_std <- NULL
  } else {
    niter <- niter*nvals
    alpha_std <- matrix(pars$alpha, ncol=1)
    pcure_std <- if (!is.null(pars$pcure)) matrix(as.vector(pars$pcure), ncol=1)
    coefs_std <- array(pars$coefs, dim=c(niter, 1, nvars))
  }

  list(alpha = alpha_std,
       coefs = coefs_std,
       pcure = pcure_std,
       niter = niter)
}
