## TODO
## use capital X for covariates not x
## Document / error that variables should be in data frame not in working env
## should loghaz be called logscale as it's not the log hazard/? but its more like a rate than a scale so that'd also be confusing.
## hazscale?  loghscale? logeta?

## prior_smooth, prior_sdnp inconsistently named

#' Flexible Bayesian parametric survival models
#'
#' Flexible Bayesian parametric survival models.  Individual data are
#' represented using M-splines and a proportional hazards or flexible
#' non-proportional hazards model.  Optionally a mixture cure version
#' of this model can be fitted.  Extrapolations can be enhanced by
#' including external aggregate data, or by including a fixed
#' background hazard in an additive hazards (relative survival) model.
#'
#' @param formula  A survival formula in standard R formula syntax, with a call to `Surv()`
#' on the left hand side.
#'
#' Covariates included on the right hand side of the formula with be modelled with
#' proportional hazards.
#'
#' @param data Data frame containing variables in `formula`.
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
#' It is intended to add facilities to the package to produce this format of data from
#' other common forms of external data, e.g. registry or population data, or elicited judgements
#' about conditional survival
#'
#' @param smooth_sd Smoothing parameter estimation.
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
#' @param cure If `TRUE`, a mixture cure model is used, where the "uncured" survival is defined by the
#' M-spline model, and the cure probability is estimated.
#'
#' @param nonprop If \code{TRUE} then a non-proportional hazards model is fitted.
#' This is achieved by modelling the spline basis weights in terms of the covariates.  See
#' the "methods" vignette for more details.   In models with multiple covariates, currently
#' there is no way to assume that some covariates have proportional hazards but others don't -
#' it is all or none.
#'
#' @param prior_loghaz Prior for the baseline log hazard scale parameter (`log(eta)`).
#'   This should be a call to a prior constructor function, such as
#'   `p_normal(0,1)` or `p_t(0,2,2)`.   Supported prior distribution families
#'   are normal (parameters mean and SD) and t distributions (parameters
#' location, scale and degrees of freedom).  The default is a normal distribution with
#' mean 0 and standard deviation 20.
#'
#' Note that `eta` is not in itself a hazard, but it is proportional to the hazard (see the vignette for the
#' full model specification).
#'
#' "Baseline" is defined
#'   by the continuous covariates taking a value of zero and factor covariates taking their
#'   reference level.  To use a different baseline, the data should be transformed
#'   appropriately beforehand, so that a value of zero has a different meaning.
#'
#' For continuous covariates, it helps for both computation and interpretation to define the value of zero to
#' denote a typical value in the data, e.g. the mean.
#'
#' @param prior_loghr Priors for log hazard ratios.  This should be a call to
#'   `p_normal()` or `p_t()`.  A list of calls can also be provided, to give
#'   different priors to different coefficients, where the name
#'   of each list component matches the name of the coefficient, e.g.
#'   ```list("age45-59" = p_normal(0,1), "age60+" = p_t(0,2,3))```
#'
#'   The default is `p_normal(0,2.5)` for all coefficients.
#'
#' @param prior_smooth Gamma prior for the smoothing standard deviation, in
#'   models where this is estimated.  This should be a call to `p_gamma()`.  The
#'   default is `p_gamma(2,1)`.
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
#' @param prior_sdnp Prior for the standard deviation parameters that smooth the non-proportionality
#' effects over time in non-proportional hazards models.  This should be a call to `p_gamma()`
#' or a list of calls to `p_gamma()` with one component per covariate, as in `prior_loghr`.
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
#' the "background" hazard at all times as a piecewise-constant (step) function.
#' Each row gives the "background" hazard between the specified time and the
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
#'   `iknots`: Internal knots.  If this is not supplied, then the number of
#'   knots is taken from `df`, and their location is taken from equally-spaced
#'   quantiles of the observed event times (concatenated with the distinct
#'   follow-up times in the external data).
#'
#'   `bknots`: Boundary knots.  If this is not supplied, the boundary knots are
#'   set to zero and the maximum event time (including the follow-up times in
#'   the external data).
#'
#'   `df`: Degrees of freedom, i.e. the number of parameters (or basis terms)
#'   that define the hazard as a spline function of time.  Defaults to 10.  Note
#'   this does not necessarily overfit because the function is smoothed through
#'   the prior.
#'
#' `degree`: Polynomial degree used for the basis function. The default is 3, giving a cubic.
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
##' @param loo Compute leave-one-out cross-validation statistics.  This is done by default. Set to
##' \code{FALSE} to not compute them.
##' If these statistics are computed, then they are returned in the \code{loo} component of the
##' object returned by \code{survextrap}.  See the \code{"examples"} vignette for some explanation
##' of these.
##'
#' @param ... Additional arguments to supply to control the Stan fit, passed to the appropriate
#' \pkg{rstan} function, depending on which is chosen through the `fit_method` argument.
#'
#' @return A list of objects defining the fitted model.  A particularly important component is
#'
#' `stanfit`: The object returned by `rstan` containing samples from the fitted model.  See the [rstan](https://mc-stan.org/rstan/index.html) documentation.
#'
#' @export
survextrap <- function(formula,
                       data,
                       external=NULL,
                       smooth_sd = "bayes",
                       coefs_mean = NULL,
                       cure = FALSE,
                       nonprop = FALSE,
                       prior_loghaz = p_normal(0, 20),
                       prior_loghr = NULL,
                       prior_smooth = p_gamma(2,1),
                       prior_cure = p_beta(1,1),
                       prior_logor_cure = NULL,
                       prior_sdnp = p_gamma(2,1),
                       backhaz = NULL,
                       mspline = NULL,
                       fit_method = "mcmc",
                       loo = (fit_method=="mcmc"),
                       ...)
{
    mf <- stats::model.frame(terms(formula), data)

    t_beg <- make_t(mf, type = "beg") # entry time
    t_end <- make_t(mf, type = "end") # exit  time
    t_upp <- make_t(mf, type = "upp") # upper time for interval censoring [ not implemented yet ]
    status <- make_d(mf)
    t_event <- aa(t_end[status == 1]) # exact event time
    t_rcens <- aa(t_end[status == 0]) # right censoring time
    nevent <- sum(status == 1)
    nrcens <- sum(status == 0)
    ind_event <- which(status==1)
    ind_rcens <- which(status==0)

    x <- make_x(formula,data)
    ncovs <- x$ncovs
    x_event <- x$X[ind_event, , drop = FALSE]
    x_rcens <- x$X[ind_rcens, , drop = FALSE]

    xcure <- make_xcure(cure,data=data)
    cure <- xcure$cure
    cure_formula <- xcure$cure_formula
    ncurecovs <- xcure$ncovs
    xcure_event <- xcure$X[ind_event, , drop = FALSE]
    xcure_rcens <- xcure$X[ind_rcens, , drop = FALSE]

    t_tmp <- sum(rowMeans(cbind(t_end, t_upp), na.rm = TRUE) - t_beg)
    d_tmp <- sum(!status == 0)

    if (is.data.frame(backhaz)) {
      validate_backhaz_df(backhaz)
      backhaz_event <- backhaz$hazard[findInterval(t_event, backhaz$time)]
      backhaz_df <- backhaz
    } else if (!is.null(backhaz)) {
      if (!is.character(backhaz)) stop("`backhaz` should be a data frame or a character string giving the name of a column in the data")
      bh <- data[[backhaz]]
      backhaz_event <- bh[ind_event]
      backhaz_df <- NULL
    } else backhaz_event <- backhaz_df <- NULL
    relative <- !is.null(backhaz_event)
    backhaz_event <- if (relative) aa(backhaz_event) else aa(numeric(nevent))

    external <- parse_external(external, formula, x, xcure, cure_formula, relative)
    t_ext_stop <- aa(external$stop)
    t_ext_start <- aa(external$start)
    r_ext <- aa(external$r)
    n_ext <- aa(external$n)
    nextern <- external$nextern
    x_ext <- external$X
    xcure_ext <- external$Xcure
    tmax <- max(c(t_end,t_upp,external$tmax), na.rm = TRUE)

    mspline <- make_basehaz(mspline    = mspline,
                            times          = t_end,
                            times_ext      = unique(c(t_ext_start, t_ext_stop)),
                            status         = status,
                            tmin          = min(t_beg),
                            tmax          = tmax)

    basis_event  <- make_basis(t_event, mspline)
    ibasis_event <- make_basis(t_event, mspline, integrate = TRUE)
    ibasis_rcens <- make_basis(t_rcens, mspline, integrate = TRUE)
    nvars <- mspline$nvars

    if (is.null(coefs_mean)){
        coefs_mean <- mspline_constant_weights(mspline)
    } else {
      coefs_mean <- validate_coefs_mean(coefs_mean)
    }
    b_mean <- aa(log(coefs_mean[-1] / coefs_mean[1]))

    est_smooth <- (smooth_sd == "bayes")
    if (est_smooth) smooth_sd_init <- 1

    ibasis_ext_stop <- if (nextern>0) make_basis(t_ext_stop, mspline, integrate = TRUE) else matrix(nrow=0, ncol=nvars)
    ibasis_ext_start <- if (nextern>0) make_basis(t_ext_start, mspline, integrate = TRUE) else matrix(nrow=0, ncol=nvars)

    if (is.data.frame(backhaz)) {
      backsurv_ext_stop <- exp(-get_cum_backhaz(t_ext_stop, backhaz))
      backsurv_ext_start <- exp(-get_cum_backhaz(t_ext_start, backhaz))
    } else {
      backsurv_ext_stop <- aa(external$backsurv_stop)
      backsurv_ext_start <- aa(external$backsurv_start)
    }
    if (nonprop) {
      if (ncovs==0) {
        warning("Ignoring non-proportional hazards model specification, since no covariates in model. ")
        nonprop <- FALSE
      }
    }

    priors <- get_priors(prior_loghaz, prior_loghr, prior_smooth, prior_cure, prior_logor_cure,
                         x, xcure, est_smooth, nonprop, prior_sdnp)

    standata <- nlist(nevent, nrcens, nvars, nextern, ncovs,
                      basis_event, ibasis_event, ibasis_rcens,
                      ibasis_ext_stop, ibasis_ext_start,
                      x_event, x_rcens,
                      r_ext, n_ext, x_ext,
                      ncurecovs, xcure_event, xcure_rcens, xcure_ext,
                      b_mean,
                      est_smooth,
                      cure,
                      relative, backhaz_df,
                      backsurv_ext_stop, backsurv_ext_start,
                      prior_loghaz_dist = priors$loghaz$distid,
                      prior_loghaz = as.numeric(unlist(priors$loghaz[c("location","scale","df")])),
                      prior_loghr_dist = aa(priors$loghr$distid),
                      prior_loghr_location = aa(priors$loghr$location),
                      prior_loghr_scale = aa(priors$loghr$scale),
                      prior_loghr_df = aa(priors$loghr$df),
                      prior_logor_cure_dist = aa(priors$logor_cure$distid),
                      prior_logor_cure_location = aa(priors$logor_cure$location),
                      prior_logor_cure_scale = aa(priors$logor_cure$scale),
                      prior_logor_cure_df = aa(priors$logor_cure$df),
                      prior_smooth = as.numeric(unlist(priors$smooth[c("shape","rate")])),
                      prior_cure = as.numeric(unlist(priors$cure[c("shape1","shape2")])),
                      modelid = 1,
                      nonprop,
                      prior_sdnp = priors$sdnp
                      )
    pcure_init <- if (cure) 0.5 else numeric()
    staninit <- list(gamma = aa(0),
                     loghr = aa(rep(0, standata$ncovs)),
                     beta_err = aa(rep(0, standata$nvars-1)),
                     smooth_sd = aa(if(standata$est_smooth) smooth_sd_init else numeric()),
                     pcure = aa(pcure_init))
    if (identical(smooth_sd, "eb")){
        smooth_sd <- eb_smoothness(standata, staninit, prior_smooth)
    }
    standata$smooth_sd_fixed <- if (est_smooth) aa(numeric()) else aa(smooth_sd)

    stan_optimizing_ops <- function(...){
        ops <- list(...)
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
                       fit_method, cure_formula,
                       backhaz=backhaz_df)
    standata_keep <- standata[c("nvars","ncovs","ncurecovs")]
    model_keep <- nlist(cure, est_smooth, nonprop)
    spline_keep <- nlist(mspline)
    covinfo_names <- c("xnames","xlevs","xinds","xbar","mfbase","mfzero")
    x <- list(x = x[covinfo_names])
    xcure <- list(xcure = xcure[covinfo_names])
    prior_keep <- list(priors=priors)
    prioretc_keep <- nlist(coefs_mean, smooth_sd)
    prior_pred <- get_prior_pred(mspline=mspline, 
                                 coefs_mean=coefs_mean, prior_smooth=prior_smooth,
                                 prior_loghaz=prior_loghaz, prior_loghr=prior_loghr,
                                 nonprop=nonprop, prior_sdnp=prior_sdnp)
    res <- c(misc_keep, standata_keep, model_keep, spline_keep, x, xcure,
             prior_keep, prioretc_keep, nlist(prior_pred), nlist(km))

    class(res) <- "survextrap"
    if (loo) res$loo <- loo_survextrap(res, standata)
    res
}

## Construct functions to draw from the prior predictive distributions
## for interesting quantities in a model specified with survextrap()

get_prior_pred <- function(mspline, 
                           coefs_mean, prior_smooth, prior_loghaz, prior_loghr,
                           nonprop, prior_sdnp){
  ## Samples of basic parameters defining the spline
  sample <- function(X=NULL, nsim=100){
    prior_sample(mspline=mspline, 
                 coefs_mean=coefs_mean, prior_smooth=prior_smooth,
                 prior_loghaz=prior_loghaz, prior_loghr=prior_loghr,
                 X=X, prior_sdnp=prior_sdnp,
                 nsim=nsim)
  }
  ## Summary of the hazard if coefs_mean was set to a constant hazard
  haz_const <- function(nsim=10000, quantiles=c(0.025, 0.5, 0.975)){
    prior_haz_const(mspline=mspline, 
                    prior_loghaz = prior_loghaz, 
                    nsim=nsim, quantiles=quantiles)
  }
  ## Samples of hazard curves over time
  haz <- function(X=NULL, tmin=mspline$bknots[1], tmax=mspline$bknots[2], nsim=10){
    mspline_priorpred(iknots=mspline$iknots, bknots=mspline$bknots, degree=mspline$degree, 
                      coefs_mean=coefs_mean, prior_smooth=prior_smooth,
                      prior_loghaz=prior_loghaz, prior_loghr=prior_loghr,
                      X=X, prior_sdnp=prior_sdnp,
                      tmin=tmin, tmax=tmax, 
                      nsim=nsim)
  }

  ## SD describing variation in hazard function with time
  haz_sd <- function(X=NULL, tmin=mspline$bknots[1], tmax=mspline$bknots[2], nsim=100,
                     hq = c(0.1, 0.9), 
                     quantiles = c(0.025, 0.5, 0.975)){
    prior_haz_sd(mspline=mspline, 
                 coefs_mean=coefs_mean, prior_smooth=prior_smooth,
                 prior_loghaz=prior_loghaz, prior_loghr=prior_loghr,
                 X=X, prior_sdnp=prior_sdnp,
                 tmin=tmin, tmax=tmax, nsim=nsim, hq=hq, quantiles=quantiles)
  }
  ## SD describing variation in hazard ratio function with time
  hr_sd <- function(X, X0, tmin=mspline$bknots[1], tmax=mspline$bknots[2], nsim=100,
                    quantiles = c(0.025, 0.5, 0.975)){
    prior_hr_sd(mspline=mspline, 
                coefs_mean=coefs_mean, prior_smooth=prior_smooth,
                prior_loghaz=prior_loghaz, prior_loghr=prior_loghr,
                X=X, X0=X0, prior_sdnp=prior_sdnp,
                tmin=tmin, tmax=tmax, nsim=nsim, quantiles=quantiles)
  }
  nlist(sample, haz, haz_const, haz_sd, hr_sd)
}

## No centering is done here.  It assumes that any covariates are pre-centered
## around a meaningful value.   This corresponds to the intercept that the prior is placed on

make_x <- function(formula, data, xlevs=NULL){
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
    X <- model.matrix(formula, mf, xlev = xlevs)
    X <- drop_intercept(X)
    ncovs <- NCOL(X)
    xnames <- colnames(X)

    nlist(X, N = NROW(X),
          ncovs = NCOL(X), xnames, factors, numerics, xlevs, mfbase)
}


eb_smoothness <- function(standata, staninit, prior_smooth){
    standata$est_smooth <- 1
    standata$smooth_sd_fixed  <- aa(numeric()) # dummy
    prior <- get_prior_smooth(prior_smooth, est_smooth=TRUE)
    standata$prior_smooth <- as.numeric(unlist(prior[c("shape","rate")]))
    staninit$smooth_sd <- aa(1)
    fits <- rstan::optimizing(stanmodels$survextrap, data=standata, init=staninit, hessian=FALSE, verbose=TRUE)
    if (fits$return_code==0){
        smooth_sd <- fits$par["smooth_sd[1]"]
    } else {
        warning("Empirical Bayes estimation of smoothness parameter failed, continuing with default")
        smooth_sd <- 1
    }
    smooth_sd
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

# Check object is a Surv object with a valid type
#
# @param x A Surv object. That is, the LHS of a formula as evaluated in a
#   data frame environment.
# @param ok_types A character vector giving the valid types of Surv object.
# @return A Surv object.
validate_surv <- function(x, ok_types = c("right", "counting",
                                          "interval", "interval2")) {
  if (!inherits(x, "Surv"))
    stop2("LHS of 'formula' must be a 'Surv' object.")
  if (!attr(x, "type") %in% ok_types)
    stop2("Surv object type must be one of: ", paste(ok_types, collapse=","))
  x
}


# Return the response vector (time)
#
# @param model_frame The model frame.
# @param type The type of time variable to return:
#   "beg": the entry time for the row in the survival data,
#   "end": the exit  time for the row in the survival data,
#   "gap": the difference between entry and exit times,
#   "upp": if the row involved interval censoring, then the exit time
#          would have been the lower limit of the interval, and "upp"
#          is the upper limit of the interval.
# @return A numeric vector.
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

# Return the response vector (status indicator)
#
# @param model_frame The model frame.
# @return A numeric vector.
make_d <- function(model_frame) {

  resp <- if (survival::is.Surv(model_frame))
    model_frame else model.response(model_frame)
  surv <- attr(resp, "type")
  err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")

  switch(surv,
         "right"     = as.vector(resp[, "status"]),
         stop(err))
}

parse_external <- function(external, formula, x, xcure, cure_formula, relative=FALSE){
  if (is.null(external))
    extl <- list(nextern=0, stop=numeric(), start=numeric(),
                 r=integer(), n=integer(),
                 backsurv_stop=numeric(), backsurv_start=numeric(),
                 tmax=0)
  else {
    validate_external(external, x, xcure, relative)
    extl <- c(as.list(external), nextern=nrow(external),
              tmax=max(external$stop))
    if (!relative) {
      extl$backsurv_stop <- extl$backsurv_start <- rep(1, extl$nextern)
    }
    if (x$ncovs>0){
      form <- delete.response(terms(formula))
      X <- model.matrix(form, external, xlev = x$xlevs)
      X <- drop_intercept(X)
    }
    if (xcure$ncovs>0){
      form <- delete.response(terms(cure_formula))
      Xcure <- model.matrix(form, external, xlev = xcure$xlevs)
      Xcure <- drop_intercept(Xcure)
    }
  }
  if ((extl$nextern==0) || (x$ncovs==0))
    X <- array(dim=c(extl$nextern, x$ncovs))
  if ((extl$nextern==0) || (xcure$ncovs==0))
    Xcure <- array(dim=c(extl$nextern, xcure$ncovs))
  extl <- c(extl, nlist(X), nlist(Xcure))
  extl
}

validate_external <- function(external, x, xcure, relative=FALSE){
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
  if (relative) {
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
  mspline_basis(times, iknots = mspline$iknots, bknots=mspline$bknots,
                degree = mspline$degree, integrate = integrate)
}

aa <- function(x, ...) as.array(x,...)


make_basehaz <- function(mspline,
                         times,
                         times_ext,
                         status,
                         tmin, tmax){
    df     <- mspline$df
    iknots <- mspline$iknots
    bknots <- mspline$bknots
    degree <- mspline$degree
    bh_names <- c("df","iknots","bknots","degree")
    bad_names <- setdiff(names(mspline), bh_names)
    if (length(bad_names) > 0) {
      blist <- paste(bad_names, collapse=",")
      plural <- if (bad_names > 1) "s" else ""
      warning(sprintf("Element%s `%s` of `mspline` is unused. Ignoring.", plural, blist))
    }
    if (is.null(df))
      df <- 10L
    if (is.null(degree))
      degree <- 3L # cubic splines
    if (df < degree + 2)
      stop(sprintf("df = %s and degree = %s, but df - degree should be >= 2", df, degree))
    tt <- times[status == 1] # uncensored event times
    if (is.null(iknots) && !length(tt)) {
        warning2("No observed events found in the data. Censoring times will ",
                 "be used to evaluate default knot locations for splines.")
        tt <- times
    }
    if (is.null(bknots)){
        bknots <- c(0, tmax)
    } else {
      validate_knots(bknots, "bknots")
      bknots <- sort(bknots)  # No restriction on boundary knots being outside the data, as we can use constant hazard outside the boundary knots.
#      if (bknots[1] > tmin) stop(sprintf("lower boundary knot (%s) should be <= minimum entry time (%s)", bknots[1], tmin))
#      if (bknots[2] < tmax) stop(sprintf("upper boundary knot (%s) should be >= maximum follow up time (%s)", bknots[2], tmax)) # unnecessary should it? can't we have a constant hazard?
    }
    if (!is.null(iknots)) {
      validate_knots(iknots, "iknots")
      if (!all(iknots > bknots[1])) stop(sprintf("`iknots` should all be greater than lower boundary knot (%s)", bknots[1]))
      if (!all(iknots < bknots[2])) stop(sprintf("`iknots` should all be less than upper boundary knot (%s)", bknots[2]))
    }
    names(bknots) <- c("lower","upper")
    ttk <- unique(c(tt, times_ext))
    iknots <- get_iknots(ttk, df = df, iknots = iknots, degree = degree, intercept = TRUE)

    nvars  <- df
    knots <- c(bknots[1], iknots, bknots[2])
    nlist(nvars, iknots, bknots, degree, df, knots)
}

validate_coefs_mean <- function(coefs){
  if (!is.numeric(coefs)) stop("`coefs_mean` must be numeric")
  if (!all(coefs > 0)) stop("`coefs_mean` must all be > 0")
  coefs / sum(coefs)
}

# Return a vector with internal knots for 'x', based on evenly spaced quantiles
#
# @param x A numeric vector.
# @param df The degrees of freedom. If specified, then 'df - degree - intercept'.
#   knots are placed at evenly spaced percentiles of 'x'. If 'iknots' is
#   specified then 'df' is ignored.
# @param degree Non-negative integer. The degree for the spline basis.
# @param iknots Optional vector of internal knots.
# @return A numeric vector of internal knot locations, or NULL if there are
#   no internal knots.
get_iknots <- function(x, df = 5L, degree = 3L, iknots = NULL, intercept = FALSE) {

  # obtain number of internal knots
  if (is.null(iknots)) {
    nk <- df - degree - intercept
  } else {
    nk <- length(iknots)
  }
  # validate number of internal knots
    if (nk < 0) {
    stop("Number of internal knots cannot be negative.")
  }
  # if no internal knots then return empty vector
  if (nk == 0) {
    return(numeric(0))
  }
  # obtain default knot locations if necessary
  if (is.null(iknots)) {
    iknots <- quantile(x, probs=seq(1, nk)/(nk + 1)) # evenly spaced percentiles
  }

  return(iknots)
}

make_xcure <- function(cure, data){
    if (isTRUE(cure)) {
        cure_formula <- ~1
    } else if (inherits(cure,"formula")) {
        cure_formula <- cure
    } else if (isFALSE(cure)) {
        cure_formula <- NULL
    } else stop("`cure` must either be TRUE, FALSE or a model formula")
    if (!is.null(cure_formula)){
        Terms <- terms(cure_formula)
        xcure <- make_x(cure_formula, data)
        cure <- TRUE
    } else {
        xcure <- list(ncovs = 0,
                      X = matrix(nrow=nrow(data), ncol=0))
    }
    c(nlist(cure, cure_formula), xcure)
}
