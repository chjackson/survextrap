#' Flexible Bayesian parametric survival models
#'
#' Flexible Bayesian parametric survival models.  Individual data are
#' represented using M-splines and a proportional hazards or flexible
#' non-proportional hazards model.  External aggregate data can be
#' included, for example, to enable extrapolation outside the
#' individual data.  A fixed background hazard can also be included in
#' an additive hazards (relative survival) model.  Mixture cure
#' versions of these models can also be used.
#'
#' @param formula  A survival formula in standard R formula syntax, with a call to `Surv()`
#' on the left hand side.
#'
#' Covariates included on the right hand side of the formula with be modelled with
#' proportional hazards.
#'
#' @param data Data frame containing variables in `formula`.  Variables should be
#' in a data frame, and not in the working environment.
#'
#' @param external External data as a data frame of aggregate survival counts with columns:
#'
#' `start`: Start time
#'
#' `stop`: Follow-up time
#'
#' `n`: Number of people alive at `start`
#'
#' `r`: Number of those people who are still alive at `stop`
#'
#' @param cure If `TRUE`, a mixture cure model is used, where the "uncured" survival is defined by the
#' M-spline model, and the cure probability is estimated.
#'
#' @param nonprop If \code{TRUE} then a non-proportional hazards model is fitted.
#' This is achieved by modelling the spline basis coefficients in terms of the covariates.  See
#' the [methods vignette](https://chjackson.github.io/survextrap/articles/methods.html) for more details.   In models with multiple covariates, currently
#' there is no way to assume that some covariates have proportional hazards but others don't -
#' it is all or none.
#'
#' @param prior_hscale Prior for the baseline log hazard scale parameter (`log(eta)`).
#'   This should be a call to a prior constructor function, such as
#'   `p_normal(0,1)` or `p_t(0,2,2)`.   Supported prior distribution families
#'   are normal (parameters mean and SD) and t distributions (parameters
#' location, scale and degrees of freedom).  The default is a normal distribution with
#' mean 0 and standard deviation 20.
#'
#' Note that `eta` is not in itself a hazard, but it is proportional to the hazard (see the vignette for the
#' full model specification).
#'
#' "Baseline" is defined by the continuous covariates taking a value
#'   of zero and factor covariates taking their reference level.  To
#'   use a different baseline, the data should be transformed
#'   appropriately beforehand, so that a value of zero has a different
#'   meaning.  For continuous covariates, it helps for both
#'   computation and interpretation to define the value of zero to
#'   denote a typical value in the data, e.g. the mean.
#'
#' @param prior_loghr Priors for log hazard ratios.  This should be a call to
#'   `p_normal()` or `p_t()`.  A list of calls can also be provided, to give
#'   different priors to different coefficients, where the name
#'   of each list component matches the name of the coefficient, e.g.
#'   ```list("age45-59" = p_normal(0,1), "age60+" = p_t(0,2,3))```
#'
#'   The default is `p_normal(0,2.5)` for all coefficients.
#'
#' @param prior_hsd Gamma prior for the standard deviation that
#'   controls the variability over time (or smoothness) of the hazard
#'   function.  This should be a call to `p_gamma()`.  The default is
#'   `p_gamma(2,1)`.  See \code{\link{prior_haz_sd}} for a way to
#'   calibrate this to represent a meaningful belief.
#'
#' @param prior_cure Prior for the baseline cure probability.  This should be a
#'   call to `p_beta()`.  The default is a uniform prior, `p_beta(1,1)`.
#'   Baseline is defined by the mean of continuous covariates and the reference
#'   level of factor covariates.
#'
#' @param prior_logor_cure Priors for log odds ratios on cure probabilities.
#'   This should be a call to `p_normal()` or `p_t()`.  The default is
#'   `p_normal(0,2.5)`.
#'
#' @param prior_hrsd Prior for the standard deviation parameters that
#'   smooth the non-proportionality effects over time in
#'   non-proportional hazards models.  This should be a call to
#'   `p_gamma()` or a list of calls to `p_gamma()` with one component
#'   per covariate, as in `prior_loghr`.  See
#'   \code{\link{prior_hr_sd}} for a way to calibrate this to
#'   represent a meaningful belief.
#'
#' @param backhaz Background hazard, that is, for causes of death other than the
#'   cause of interest. This defines a "relative survival" model where the
#'   overall hazard is the sum of a cause-specific hazard and a background
#'   hazard.   The background hazard is assumed to be known, and the
#'   cause-specific hazard is modelled with the flexible parametric model.
#'
#' The background hazard can be supplied in two forms.  The meaning of predictions
#' from the model depends on which of these is used.
#'
#' (a) A data frame with columns \code{"hazard"} and \code{"time"}, specifying
#' the background hazard at all times as a piecewise-constant (step) function.
#' Each row gives the background hazard between the specified time and the
#' next time.  The first element of \code{"time"} should be 0, and the final row
#' specifies the hazard at all times greater than the last element of
#' \code{"time"}.
#'
#' Predictions from the modelled fitted by `survextrap` will then include this
#' background hazard, because it is known at all times.
#'
#' (b) The (quoted) name of a variable in the data giving the background
#' hazard.  For censored cases, the exact value does
#' not matter.  The predictions from `survextrap` will
#' then describe the excess hazard or survival on top of this background.
#' The overall hazard cannot be predicted in general, because the background hazard is
#' only specified over a limited range of time.
#'
#' If there is external data and `backhaz` is supplied in form (b), then the
#' user should also supply the background survival at the start and stop points
#' in columns of the external data named `"backsurv_start"` and
#' `"backsurv_stop"`.  This should describe the same reference population as
#' `backhaz`, though the package does not check for consistency between these.
#'
#' Leave-one-out cross-validation currently does not take into account any
#' background hazard (effectively assuming it to be zero).
#'
#' @param mspline A list of control parameters defining the spline model.
#'
#'   `knots`: Spline knots.  If this is not supplied, then the number
#'   of knots is taken from `df`, and their location is taken from
#'   equally-spaced quantiles of the observed event times in the
#'   individual-level data.
#'
#'   `add_knots`: Any extra knots beyond those chosen from the
#'   individual-level data (or supplied in `knots`).  This is intended
#'   to be used when there are \code{external} data included in the
#'   model.  External data are typically outside the time period
#'   covered by the individual data.  `add_knots` would then be chosen
#'   to span the time period covered by the external data, so that the
#'   hazard trajectory can vary over that time.
#'
#'   If there are external data, and both `knots` and `add_knots` are
#'   omitted, then a default set of knots is chosen to span both the
#'   individual and external data, by taking the quantiles of a vector
#'   defined by concatenating the individual-level event times with
#'   the \code{start} and \code{stop} times in the external data.
#'
#'   `df`: Degrees of freedom, i.e. the number of parameters (or basis
#'    terms) intended to result from choosing knots based on quantiles
#'    of the data.  The total number of parameters will then be `df`
#'    plus the number of additional knots specified in
#'    `add_knots`. `df` defaults to 10.  This does not necessarily
#'    overfit, because the function is smoothed through the prior.
#'
#'   `degree`: Polynomial degree used for the basis function. The
#'   default is 3, giving a cubic.
#'
#'   `bsmooth`: If \code{TRUE} (on by default) the spline is smoother
#'    at the highest knot, by defining the derivative and second derivative
#'    at this point to be zero.
#'
#' @param hsd Smoothing parameter estimation.
#'
#' `"bayes"`: the smoothing parameter is estimated by full Bayes.
#'
#'  `"eb"`: empirical Bayes is used.
#'
#' Alternatively, if a number is supplied here, then the smoothing parameter is fixed to this number.
#'
#' @param coefs_mean Spline basis coefficients that define the prior mean for the hazard function. By
#' default, these are set to values that define a constant hazard function.  They are normalised to
#' sum to 1 internally (if they do not already).
#'
#' @param fit_method Method from \pkg{rstan} used to fit the model.  Defaults to MCMC.
#'
#'  If \code{fit_method="mcmc"} then a sample from the posterior is drawn using Markov Chain Monte Carlo
##' sampling, via \pkg{rstan}'s \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}} function.
##' This is the default.  It is the most accurate but the slowest method.
##'
##'   If \code{fit_method="opt"}, then instead of an MCMC sample from the posterior,
##'   `survextrap` returns the posterior mode calculated using optimisation, via
##'   \pkg{rstan}'s \code{\link[rstan:stanmodel-method-optimizing]{rstan::optimizing()}} function.
##'   A sample from a normal approximation to the (real-line-transformed)
##'   posterior distribution is drawn in order to obtain credible intervals.
##'
##'   If \code{fit_method="vb"}, then variational Bayes methods are used, via \pkg{rstan}'s
##'   \code{\link[rstan:stanmodel-method-vb]{rstan::vb()}} function.  This is labelled as "experimental" by
##'   \pkg{rstan}.  It might give a better approximation to the posterior
##'   than \code{fit_method="opt"}, but has not been investigated in depth for these models.
##'
##' @param loo Compute leave-one-out cross-validation statistics.
##'   This is done by default. Set to \code{FALSE} to not compute
##'   them.  If these statistics are computed, then they are returned
##'   in the \code{loo} and \code{loo_external} components of the
##'   object returned by \code{survextrap}.  \code{loo} describes the
##'   fit of the model to the individual-level data, and
##'   \code{loo_external} describes fit to the external data.
##' 
##' See the \code{"examples"} vignette for more explanation of these.
##'
#' @param ... Additional arguments to supply to control the Stan fit, passed to the appropriate
#' \pkg{rstan} function, depending on which is chosen through the `fit_method` argument.
#'
#' @return A list of objects defining the fitted model.  These are not
#'   intended to be extracted directly by users.  Instead see
#'   \code{\link{summary.survextrap}} for summarising the parameter
#'   estimates, and the functions \code{\link{hazard}},
#'   \code{\link{survival}}, \code{\link{rmst}} and others for
#'   computing interesting summaries of the fitted survival
#'   distribution.
#'
#' The object returned by `rstan` containing samples from the fitted
#' model is returned in the `stanfit` component.  See the
#' [rstan](https://mc-stan.org/rstan/index.html) documentation.  The
#' function \code{\link{get_draws}} is provided to convert this to a
#' simple matrix of posterior samples with all chains collapsed.
#'
#' @export
survextrap <- function(formula,
                       data,
                       external=NULL,
                       cure = FALSE,
                       nonprop = FALSE,
                       prior_hscale = p_normal(0, 20),
                       prior_loghr = NULL,
                       prior_hsd = p_gamma(2,1),
                       prior_cure = p_beta(1,1),
                       prior_logor_cure = NULL,
                       prior_hrsd = p_gamma(2,1),
                       backhaz = NULL,
                       mspline = NULL,
                       hsd = "bayes",
                       coefs_mean = NULL,
                       fit_method = "mcmc",
                       loo = (fit_method=="mcmc"),
                       ...)
{
    mf <- stats::model.frame(terms(formula), data)
    td <- make_td(mf)
    x <- make_x(formula, data, td)
    xcure <- make_xcure(cure, data, td)
    backhaz <- make_backhaz(backhaz, data, td)
    external <- make_external(external, formula, x, xcure, backhaz)
    mspline <- make_mspline(mspline, td, external)

    basis_event  <- make_basis(td$t_event, mspline)
    ibasis_event <- make_basis(td$t_event, mspline, integrate = TRUE)
    ibasis_rcens <- make_basis(td$t_rcens, mspline, integrate = TRUE)
    nvars <- mspline$nvars

    if (is.null(coefs_mean)){
        coefs_mean <- mspline_constant_coefs(mspline)
    } else {
      coefs_mean <- validate_coefs_mean(coefs_mean)
    }
    b_mean <- aa(log(coefs_mean[-1] / coefs_mean[1]))

    est_hsd <- (hsd == "bayes")
    if (est_hsd) hsd_init <- 1

    ibasis_ext_stop <- if (external$nextern>0) make_basis(external$stop, mspline, integrate = TRUE) else matrix(nrow=0, ncol=nvars)
    ibasis_ext_start <- if (external$nextern>0) make_basis(external$start, mspline, integrate = TRUE) else matrix(nrow=0, ncol=nvars)

    if (nonprop) {
      if (x$ncovs==0) {
        warning("Ignoring non-proportional hazards model specification, since no covariates in model. ")
        nonprop <- FALSE
      }
    }

    priors <- get_priors(prior_hscale, prior_loghr, prior_hsd, prior_cure, prior_logor_cure,
                         x, xcure, est_hsd, nonprop, prior_hrsd)

    standata <- nlist(nevent=td$nevent, nrcens=td$nrcens,
                      nvars, nextern=external$nextern, ncovs = x$ncovs,
                      basis_event, ibasis_event, ibasis_rcens,
                      ibasis_ext_stop, ibasis_ext_start,
                      x_event=x$event, x_rcens=x$rcens,
                      r_ext = external$r,
                      n_ext = external$n,
                      x_ext = external$X,
                      ncurecovs = xcure$ncovs,
                      xcure_event = xcure$event, xcure_rcens = xcure$rcens,
                      xcure_ext = external$Xcure,
                      b_mean,
                      est_hsd,
                      cure = xcure$cure,
                      relative=backhaz$relative,
                      backhaz_event = backhaz$event,
                      backsurv_ext_stop = external$backsurv_stop,
                      backsurv_ext_start = external$backsurv_start,
                      prior_hscale_dist = priors$hscale$distid,  # change name in stan too?
                      prior_hscale = as.numeric(unlist(priors$hscale[c("location","scale","df")])),
                      prior_loghr_dist = aa(priors$loghr$distid),
                      prior_loghr_location = aa(priors$loghr$location),
                      prior_loghr_scale = aa(priors$loghr$scale),
                      prior_loghr_df = aa(priors$loghr$df),
                      prior_logor_cure_dist = aa(priors$logor_cure$distid),
                      prior_logor_cure_location = aa(priors$logor_cure$location),
                      prior_logor_cure_scale = aa(priors$logor_cure$scale),
                      prior_logor_cure_df = aa(priors$logor_cure$df),
                      prior_hsd = as.numeric(unlist(priors$hsd[c("shape","rate")])),
                      prior_cure = as.numeric(unlist(priors$cure[c("shape1","shape2")])),
                      modelid = 1,
                      nonprop,
                      prior_hrsd = priors$sdnp
                      )
    pcure_init <- if (xcure$cure) 0.5 else numeric()
    staninit <- list(gamma = aa(0),
                     loghr = aa(rep(0, standata$ncovs)),
                     beta_err = aa(rep(0, standata$nvars-1)),
                     hsd = aa(if(standata$est_hsd) hsd_init else numeric()),
                     pcure = aa(pcure_init))
    if (identical(hsd, "eb")){
        hsd <- eb_smoothness(standata, staninit, prior_hsd)
    }
    standata$hsd_fixed <- if (est_hsd) aa(numeric()) else aa(hsd)

    stan_optimizing_ops <- function(...){
        ops <- list(...)
        ops <- ops[names(ops) %in% .rstan_optimizing_args]
        if (is.null(ops$hessian)) ops$hessian <- TRUE
        if (is.null(ops$draws)) ops$draws <- 2000
        ops
    }
    stan_sampling_ops <- function(...){
        ops <- list(...)
        ops
    }
    stan_vb_ops <- function(...){
      ops <- list(...)
      ops <- ops[names(ops) %in% .rstan_vb_args]
      ops
    }

    stanmod <- "survextrap"
    if (fit_method=="opt")
        fits <- do.call(rstan::optimizing,
                        c(list(object=stanmodels[[stanmod]], data=standata, init=staninit),
                          stan_optimizing_ops(...)))
    else if (fit_method=="mcmc")
        fits <- do.call(rstan::sampling,
                        c(list(object=stanmodels[[stanmod]], data=standata,
                               pars = "beta", include=FALSE), stan_sampling_ops(...)))
    else if (fit_method=="vb")
        fits <- do.call(rstan::vb,
                        c(list(object=stanmodels[[stanmod]], data=standata,
                               pars = "beta", include=FALSE), stan_vb_ops(...)))
    else stop(sprintf("Unknown fit_method: %s",fit_method))

    km <- survminer::surv_summary(survfit(formula, data=data), data=data)

    misc_keep <- nlist(formula, stanfit=fits,
                       fit_method,
                       cure_formula = xcure$cure_formula,
                       backhaz=backhaz$df)
    standata_keep <- standata[c("nvars","ncovs","ncurecovs","nevent","nrcens","nextern")]
    model_keep <- nlist(cure=xcure$cure, est_hsd, nonprop)
    spline_keep <- nlist(mspline)
    covinfo_names <- c("xnames","xlevs","xinds","mfbase")
    x <- list(x = x[covinfo_names])
    xcure <- list(xcure = xcure[covinfo_names])
    prior_keep <- list(priors=priors)
    prioretc_keep <- nlist(coefs_mean, hsd)
    prior_pred <- get_prior_pred(mspline=mspline,
                                 coefs_mean=coefs_mean, prior_hsd=prior_hsd,
                                 prior_hscale=prior_hscale, prior_loghr=prior_loghr,
                                 nonprop=nonprop, prior_hrsd=prior_hrsd)
    res <- c(misc_keep, standata_keep, model_keep, spline_keep, x, xcure,
             prior_keep, prioretc_keep, nlist(prior_pred), nlist(km))

    class(res) <- "survextrap"
    if (loo) {
      res$loo <- loo_survextrap(res, standata, loglik_ipd)
      if (external$nextern > 0)
        res$loo_external <- loo_survextrap(res, standata, loglik_external)
    }
    else {
      res$loo <- sprintf("Cross-validation statistics not available with `fit_method=\"%s\"`", fit_method)
      class(res$loo) <- "message"
    }
    res
}

##' @export
print.message <- function(x, ...){
  message(x)
}

##' Make default M-spline knot specification.
##'
##' Choose default M-spline knot locations given a dataset and desired
##' number of spline parameters.  Assumes a cubic spline, and knots
##' based on quantiles of event times observed in the individual data.
##'
##' If there are also external data, then these are based on quantiles
##' of a vector defined by concatenating the event times in the
##' individual data with the unique start and stop times in the
##' external data.
##'
##' @inheritParams survextrap
##'
##' @inheritParams mspline_args
##'
##' @param add_knots Additional knots, other than those determined from the quantiles of the individual data.
##' Typically used to add a maximum knot at the time that we want to extrapolate to.
##'
##' @return A list with components
##'
##' \code{knots} Knot locations.  The number of
##' knots will be equal to \code{df} + \code{degree} + 2.
##' \code{degree} Spline polynomial degree (i.e. 3)
##' \code{nvars} Number of basis variables (an alias for \code{df})
##'
##' @export
mspline_spec <- function(formula, data, cure=FALSE, backhaz=NULL, external=NULL,
                         df=10, add_knots=NULL, degree=3, bsmooth=TRUE){
  mf <- stats::model.frame(terms(formula), data)
  td <- make_td(mf)
  x <- make_x(formula, data, td)
  xcure <- make_xcure(cure, data, td)
  backhaz <- make_backhaz(backhaz, data, td)
  external <- make_external(external, formula, x, xcure, backhaz)
  mspline <- list(df=df,degree=degree,bsmooth=bsmooth,add_knots=add_knots)
  mspline <- make_mspline(mspline, td, external)
  mspline
}

## Construct functions to draw from the prior predictive distributions
## for interesting quantities in a model specified with survextrap()

get_prior_pred <- function(mspline,
                           coefs_mean, prior_hsd, prior_hscale, prior_loghr,
                           nonprop, prior_hrsd){
  ## Samples of basic parameters defining the spline
  sample <- function(X=NULL, nsim=100){
    prior_sample(mspline=mspline,
                 coefs_mean=coefs_mean, prior_hsd=prior_hsd,
                 prior_hscale=prior_hscale, prior_loghr=prior_loghr,
                 X=X, prior_hrsd=prior_hrsd,
                 nsim=nsim)
  }
  ## Summary of the hazard if coefs_mean was set to a constant hazard
  haz_const <- function(nsim=10000, quantiles=c(0.025, 0.5, 0.975)){
    prior_haz_const(mspline=mspline,
                    prior_hscale = prior_hscale,
                    nsim=nsim, quantiles=quantiles)
  }
  ## Samples of hazard curves over time
  haz <- function(X=NULL, tmin=0, tmax=max(mspline$knots), nsim=10){
    mspline_priorpred(knots=mspline$knots, degree=mspline$degree,
                      coefs_mean=coefs_mean, prior_hsd=prior_hsd,
                      prior_hscale=prior_hscale, prior_loghr=prior_loghr,
                      X=X, prior_hrsd=prior_hrsd,
                      tmin=tmin, tmax=tmax,
                      nsim=nsim)
  }

  ## SD describing variation in hazard function with time
  haz_sd <- function(X=NULL, tmin=0, tmax=max(mspline$knots), nsim=100,
                     hq = c(0.1, 0.9),
                     quantiles = c(0.025, 0.5, 0.975)){
    prior_haz_sd(mspline=mspline,
                 coefs_mean=coefs_mean, prior_hsd=prior_hsd,
                 prior_hscale=prior_hscale, prior_loghr=prior_loghr,
                 X=X, prior_hrsd=prior_hrsd,
                 tmin=tmin, tmax=tmax, nsim=nsim, hq=hq, quantiles=quantiles)
  }
  ## SD describing variation in hazard ratio function with time
  hr_sd <- function(X, X0, tmin=0, tmax=max(mspline$knots), nsim=100,
                    quantiles = c(0.025, 0.5, 0.975)){
    prior_hr_sd(mspline=mspline,
                coefs_mean=coefs_mean, prior_hsd=prior_hsd,
                prior_hscale=prior_hscale, prior_loghr=prior_loghr,
                X=X, X0=X0, prior_hrsd=prior_hrsd,
                tmin=tmin, tmax=tmax, nsim=nsim, quantiles=quantiles)
  }
  nlist(sample, haz, haz_const, haz_sd, hr_sd)
}

## No centering is done here.  It assumes that any covariates are pre-centered
## around a meaningful value.   This corresponds to the intercept that the prior is placed on

make_x <- function(formula, data, td, xlevs=NULL){
    ## model frame containing x not factor(x), log(x) etc
    mforig <- get_all_vars(delete.response(terms(formula)), data)
    if (ncol(mforig) > 0){
        factors <- names(mforig)[sapply(mforig, is.factor)]
        numerics <- names(mforig)[sapply(mforig, is.numeric)]
        tofactors <- names(mforig)[sapply(mforig, function(x){is.character(x)|is.logical(x)})]
        for (i in tofactors)
          mforig[[i]] <- factor(mforig[[i]])
        factors <- c(factors, tofactors)
        xlevs <- lapply(mforig[factors], levels)
        if ((length(factors)==1) && (length(numerics)==0)){
            ref_levs <- lapply(xlevs, function(x)x[1])
            base_levs <- xlevs
        }
        else if (length(factors) > 0){
            ref_levs <- lapply(xlevs, function(x)x[1])
            base_levs <- ref_levs
        } else base_levs <- ref_levs <- NULL
        numeric_zeros <- lapply(mforig[numerics], function(x)0)
        mfbase <- as.data.frame(c(numeric_zeros, base_levs))
    } else factors <- numerics <- xlevs <- mfbase <- NULL

    ## model frame containing  factor(x), log(x) etc not x
    mf <- stats::model.frame(formula, data)
    xlevs <- .getXlevels(terms(mf), mf)
    X <- model.matrix(formula, mf, xlev = xlevs)
    X <- drop_intercept(X)
    ncovs <- NCOL(X)
    xnames <- colnames(X)

    event <- X[td$ind_event, , drop = FALSE]
    rcens <- X[td$ind_rcens, , drop = FALSE]

    nlist(X, N = NROW(X),
          ncovs = NCOL(X), xnames, factors, numerics, xlevs,
          ## mfbase is used as the default "newdata" in output functions
          ## i.e. zero values of continuous covariates, first levels of factors
          mfbase,
          event, rcens)
}


eb_smoothness <- function(standata, staninit, prior_hsd){
    standata$est_hsd <- 1
    standata$hsd_fixed  <- aa(numeric()) # dummy
    prior <- get_prior_hsd(prior_hsd, est_hsd=TRUE)
    standata$prior_hsd <- as.numeric(unlist(prior[c("shape","rate")]))
    staninit$hsd <- aa(1)
    fits <- rstan::optimizing(stanmodels$survextrap, data=standata, init=staninit, hessian=FALSE, verbose=TRUE)
    if (fits$return_code==0){
        hsd <- fits$par["hsd[1]"]
    } else {
        warning("Empirical Bayes estimation of smoothness parameter failed, continuing with default")
        hsd <- 1
    }
    hsd
}

validate_formula <- function(formula, needs_response = TRUE) {

  if (!inherits(formula, "formula")) {
    stop2("'formula' must be a formula.")
  }

  if (needs_response) {
    len <- length(formula)
    if (len < 3) {
      stop2("'formula' must contain a response.")
    }
  }
  as.formula(formula)
}


make_td <- function(mf){
  resp <- if (survival::is.Surv(mf)) mf else model.response(mf)
  surv <- attr(resp, "type")
  err  <- paste0("Cannot handle '", surv, "' Surv objects.")
  t_end <- switch(surv,  "right" = as.vector(resp[, "time"]),  stop(err))
  if (any(t_end<0)) stop("Some survival times are negative")
  status <- make_d(mf)
  t_event <- aa(t_end[status == 1]) # exact event time
  t_rcens <- aa(t_end[status == 0]) # right censoring time
  nevent <- sum(status == 1)
  nrcens <- sum(status == 0)
  ind_event <- which(status==1)
  ind_rcens <- which(status==0)
  nlist(t_end, status,
        t_event, t_rcens, nevent, nrcens,
        ind_event, ind_rcens)
}

make_t <- function(model_frame, type = c("beg", "end", "gap", "upp")) {
    type <- match.arg(type)
    resp <- if (survival::is.Surv(model_frame))
              model_frame else model.response(model_frame)
    surv <- attr(resp, "type")
    err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")

    t_beg <- switch(surv,
                    "right"     = rep(0, nrow(model_frame)),
                    stop(err))

    t_end <- switch(surv,
                    "right"     = as.vector(resp[, "time"]),
                    stop(err))
    if (any(t_end<0)) stop("Some survival times are negative")

    t_upp <- switch(surv,
                    "right"     = rep(NaN, nrow(model_frame)),
                    stop(err))

    switch(type,
           "beg" = t_beg,
           "end" = t_end,
           "gap" = t_end - t_beg,
           "upp" = t_upp,
           stop("Bug found: cannot handle specified 'type'."))
}

make_d <- function(model_frame) {

  resp <- if (survival::is.Surv(model_frame))
    model_frame else model.response(model_frame)
  surv <- attr(resp, "type")
  err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")

  switch(surv,
         "right"     = as.vector(resp[, "status"]),
         stop(err))
}

make_external <- function(external, formula, x, xcure, backhaz){
  if (is.null(external))
    extl <- list(nextern=0, stop=numeric(), start=numeric(),
                 r=integer(), n=integer(),
                 backsurv_stop=numeric(), backsurv_start=numeric(),
                 tmax=0)
  else {
    validate_external(external, x, xcure, backhaz)
    extl <- c(as.list(external), nextern=nrow(external),
              tmax=max(external$stop))
    if (backhaz$relative) {
      if (is.data.frame(backhaz$df)) {
        extl$backsurv_stop <- aa(exp(-get_cum_backhaz(external$stop, backhaz$df)))
        extl$backsurv_start <- aa(exp(-get_cum_backhaz(external$start, backhaz$df)))
      } else {
        extl$backsurv_stop <- aa(external$backsurv_stop)
        extl$backsurv_start <- aa(external$backsurv_start)
      }
    } else {
      extl$backsurv_stop <- extl$backsurv_start <- aa(rep(1, extl$nextern))
    }
    if (x$ncovs>0){
      form <- delete.response(terms(formula))
      mf <- model.frame(form, external, xlev = x$xlevs)
      X <- drop_intercept(model.matrix(form, mf))
    }
    if (xcure$ncovs>0){
      form <- delete.response(terms(xcure$cure_formula))
      mf <- model.frame(formula, external, xlev = xcure$xlevs)
      Xcure <- drop_intercept(model.matrix(form, mf))
    }
  }
  if ((extl$nextern==0) || (x$ncovs==0))
    X <- array(dim=c(extl$nextern, x$ncovs))
  if ((extl$nextern==0) || (xcure$ncovs==0))
    Xcure <- array(dim=c(extl$nextern, xcure$ncovs))
  extl$stop <- aa(extl$stop)
  extl$start <- aa(extl$start)
  extl$r <- aa(extl$r)
  extl$n <- aa(extl$n)
  extl <- c(extl, nlist(X), nlist(Xcure))
  extl
}

validate_external <- function(external, x, xcure, backhaz){
  reqd_cols <- c("start","stop","r","n")
  missing_cols <- setdiff(reqd_cols, names(external))
  if (length(missing_cols) > 0) {
    mstr <- paste(paste0("\"",missing_cols,"\"", collapse=","))
    plural <- if (length(missing_cols) > 1) "s" else ""
    stop(sprintf("Column%s %s not found in `external`",plural,mstr))
  }
  covnames <- union(names(x$mfbase), names(xcure$mfbase))
  missing_cols <- setdiff(covnames, names(external))
  if (length(missing_cols) > 0) {
      mstr <- paste(paste0("\"",missing_cols,"\"", collapse=","))
      plural <- if (length(missing_cols) > 1) "s" else ""
      stop(sprintf("Covariate%s %s not found in `external`",plural,mstr))
  }
  if (backhaz$relative && !is.data.frame(backhaz$df)) {
    reqd_cols <- c("backsurv_start","backsurv_stop")
    missing_cols <- setdiff(reqd_cols, names(external))
    if (length(missing_cols) > 0) {
      mstr <- paste(paste0("\"",missing_cols,"\"", collapse=","))
      plural <- if (length(missing_cols) > 1) "s" else ""
      stop(sprintf("Column%s %s not found in `external`. These are required for relative survival models.",
                   plural,mstr))
    }
  }
}

make_basis <- function(times, mspline, integrate = FALSE) {
  N <- length(times)
  K <- mspline$nvars
  if (!N) { # times is NULL or empty vector
    return(matrix(0, 0, K))
  }
  mspline_basis(times, knots = mspline$knots, degree = mspline$degree,
                integrate = integrate, bsmooth=mspline$bsmooth)
}

aa <- function(x, ...) as.array(x,...)


make_mspline <- function(mspline, td, external){
    mspline <- mspline_default(mspline)
    knots <- sort(mspline$knots)
    df     <- mspline$df
    degree <- mspline$degree
    bsmooth <- mspline$bsmooth
    tmax <- max(c(td$t_end,external$tmax), na.rm = TRUE)
    if (!is.null(knots)) {
      validate_knots(knots, "knots")
      df_correct <- if (bsmooth) length(knots) + degree - 2 else length(knots) + degree
      if (!is.null(df)) {
        if (df != df_correct)
          stop(sprintf("df should equal %s given knots and bsmooth=%s",df_correct,bsmooth))
      } else df <- df_correct
    }
    else if (is.null(df))
      df <- 10L
    else
      if (df < degree + 2){
        ## unsmooth: order k = d+1, degree d = k-1,  df = nik+k = nik+d+1 = nk+d
        ## smooth: df - degree = nik - 1 unsmooth: df - degree = nik + 1
        minddg <- if (bsmooth) 0 else 2
        stop(sprintf("df = %s and degree = %s, but df - degree should be >= %s", df, degree, minddg))
      }
    tt <- td$t_end[td$status == 1] # uncensored event times
    if (is.null(knots) && !length(tt)) {
        warning2("No observed events found in the data. Censoring times will ",
                 "be used to evaluate default knot locations for splines.")
        tt <- td$t_end
    }
    ttk <- unique(tt)
    if (!is.numeric(mspline$add_knots))
      ttk <- sort(c(ttk, unique(c(external$start, external$stop))))
    knots <- default_knots(ttk, df = df, knots = knots,
                           degree = degree, bsmooth = bsmooth)
    if (is.numeric(mspline$add_knots)) {
      add_knots <- setdiff(unique(mspline$add_knots), knots)
      knots <- sort(c(knots, add_knots))
      df <- df + length(add_knots)
    }
    nvars  <- df
    nlist(nvars, knots, degree, bsmooth, df)
}

validate_coefs_mean <- function(coefs){
  if (!is.numeric(coefs)) stop("`coefs_mean` must be numeric")
  if (!all(coefs > 0)) stop("`coefs_mean` must all be > 0")
  coefs / sum(coefs)
}

# Return a vector with knots for 'x', based on evenly spaced quantiles
#
# @param x A numeric vector.
# @param df The degrees of freedom. If specified, then 'df - degree - 1'.
#   knots are placed at evenly spaced percentiles of 'x'. If 'iknots' is
#   specified then 'df' is ignored.
# @param degree Non-negative integer. The degree for the spline basis.
# @param knots Optional vector of knots.
# @return A numeric vector of knot locations, or NULL if there are
#   no internal knots.
default_knots <- function(x, df = 5L, degree = 3L, knots = NULL, bsmooth = TRUE) {
  # obtain number of internal knots
  if (is.null(knots)) {
    nk <- if (bsmooth) df - degree + 2 else df - degree
  } else {
    nk <- length(knots)
  }
  # validate number of knots
    if (nk < 0) {
    stop("Number of knots cannot be negative.")
  }
  # if no knots then return empty vector
  if (nk == 0) {
    return(numeric(0))
  }
  # obtain default knot locations if necessary
  if (is.null(knots)) {
    knots <- quantile(x, probs=seq(1, nk)/nk) # evenly spaced percentiles
  }
  return(knots)
}

make_xcure <- function(cure, data, td){
    if (isTRUE(cure)) {
        cure_formula <- ~1
    } else if (inherits(cure,"formula")) {
        cure_formula <- cure
    } else if (isFALSE(cure)) {
        cure_formula <- NULL
    } else stop("`cure` must either be TRUE, FALSE or a model formula")
    if (!is.null(cure_formula)){
        Terms <- terms(cure_formula)
        xcure <- make_x(cure_formula, data, td)
        cure <- TRUE
    } else {
        xcure <- list(ncovs = 0,
                      X = matrix(nrow=nrow(data), ncol=0),
                      event = matrix(nrow=td$nevent,ncol=0),
                      rcens = matrix(nrow=td$nrcens,ncol=0))
    }
    res <- c(nlist(cure, cure_formula), xcure)
    res
}



## Named arguments to Stan computation functions

## These are specified explicitly, so that the user can pass e.g.
## chains=NULL to fit_method="opt" and it will be ignored rather than
## giving an error.

## Ideally we would be able to read these from rstan by using
## formals(), but they are hidden deeply inside methods.  Instead they
## are taken from help(rstan::optimizing) - should ensure it is up to
## date

.rstan_optimizing_args <- c("object","data","seed","init","check_data","sample_file","algorithm",
                            "verbose","hessian","as_vector","draws","constrained","importance_resampling",
                            "iter","save_iterations","refresh","init_alpha","tol_obj","tol_rel_obj",
                            "tol_grad","tol_rel_grad","tol_param","history_size")

.rstan_vb_args <- c("object","data","pars","include", "seed","init","sample_file","algorithm",
                    "importance_resampling","keep_every",
                    "iter","grad_samples","elbo_samples","eta","adapt_engaged",
                    "tol_rel_obj","eval_elbo","output_samples","adapt_iter")
