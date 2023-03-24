#' Prior distributions and options
#'
#' @name priors
#'
#' @description The functions described on this page are used to specify the prior distributions for the parameters in a survextrap model.
#'
#' @param location Prior location. For the normal distribution, this is the mean.
#' Defaults to 0
#'
#' @param scale Prior scale.  For the normal distribution, this is the standard deviation.
#' Defaults to 2.5.
#'
#' @param df Prior degrees of freedom (only for Student t distribution).
#'
#' @param shape1 First shape parameter (for Beta distribution, defaults to 1).
#'
#' @param shape2 Second shape parameter (for Beta distribution, defaults to 1).
#'
#' @param shape Shape parameter (for Gamma distribution, defaults to 2).
#'
#' @param rate Rate parameter (for Gamma distribution, defaults to 1).
#'
#' @seealso \code{\link{survextrap}}.
#'
#' @return A named list to be used internally by the \pkg{survextrap} model fitting functions.
#'
#'
NULL


#' @rdname priors
#' @export
p_normal <- function(location = 0, scale = 2.5) {
  validate_positive_parameter(scale)
  res <- nlist(dist = "normal", distid=1, location, scale, df=1)
  res$r <- function(n)rnorm(n, mean=location, sd=scale)
  res$q <- function(p)qnorm(p, mean=location, sd=scale)
  class(res) <- "prior"
  res
}

#' @rdname priors
#' @export
p_t <- function(location = 0, scale = 2.5, df = 1) {
  validate_positive_parameter(scale)
  validate_positive_parameter(df)
  res <- nlist(dist = "t", distid=2, location, scale, df)
  res$r <- function(n){rt(n, df=df)*scale + location}
  res$q <- function(p){qt(p, df=df)*scale + location}
  class(res) <- "prior"
  res
}

#' @rdname priors
#' @export
p_beta <- function(shape1 = 1, shape2 = 1){
    validate_positive_parameter(shape1)
    validate_positive_parameter(shape2)
    res <- nlist(dist = "beta", shape1, shape2)
    res$r <- function(n)rbeta(n, shape1=shape1, shape2=shape2)
    res$q <- function(p)qbeta(p, shape1=shape1, shape2=shape2)
    class(res) <- "prior"
    res
}

#' @rdname priors
#' @export
p_gamma <- function(shape = 2, rate = 1){
    validate_positive_parameter(shape)
    validate_positive_parameter(rate)
    res <- nlist(dist = "gamma", shape, rate)
    res$r <- function(n)rgamma(n, shape=shape, rate=rate)
    res$q <- function(p)qgamma(p, shape=shape, rate=rate)
    class(res) <- "prior"
    res
}

#' Derive a normal prior for the hazard scale parameter based on a guess at survival times
#'
#' @inheritParams mspline_constant_weights
#'
#' @param median Best guess (prior median) for a typical survival time
#'
#' @param upper 95% credible limit for a survival time
#'
#' @export
p_meansurv <- function(median, upper, mspline=NULL){
  coefs_mean <- mspline_constant_weights(mspline)
  haz_std <-  hsurvmspline(x=mean(mspline$bknots), alpha=1, coefs=coefs_mean,
                           knots=c(mspline$iknots,mspline$bknots), degree=mspline$degree)
  haz_lower <- 1/upper
  haz_median <- 1/median
  logeta_lower <- log(haz_lower / haz_std)
  logeta_median <- log(haz_median / haz_std)
  logeta_sd <- (logeta_median - logeta_lower) / qnorm(0.975)
  p_normal(logeta_median, logeta_sd)
}


# internal ----------------------------------------------------------------

# Check for positive parameter (NULL is used in rstanarm)
#
# @param x The value to check.
# @return Either an error is thrown or \code{TRUE} is returned invisibly.
validate_positive_parameter <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x))
      stop(nm, " should be NULL or numeric", call. = FALSE)
    if (any(x <= 0))
      stop(nm, " should be positive", call. = FALSE)
  }
  invisible(TRUE)
}

get_prior_loghaz <- function(loghaz){
  if (is.null(loghaz)) loghaz <- p_normal(0, 20) # currently unused: default set in function headers
  else validate_prior(loghaz)
  loghaz
}

get_priors <- function(loghaz, loghr, smooth, cure, logor_cure, x, xcure, est_smooth,
                       nonprop, prior_sdnp){
  loghaz <- get_prior_loghaz(loghaz)
  loghr <- get_prior_coveffs(loghr, x, "loghr")
  smooth <- get_prior_smooth(smooth, est_smooth)
  validate_prior(cure)
  logor_cure <- get_prior_coveffs(logor_cure, xcure, "logor_cure")
  sdnp <- get_prior_sdnp(prior_sdnp, x, nonprop)
  nlist(loghaz, loghr, smooth, cure, logor_cure, sdnp)
}

validate_prior <- function(prior, priorname=NULL, element=NULL){
    if (is.null(priorname)) priorname <- deparse(substitute(prior))
    priorfullname <- sprintf("prior_%s", priorname)

    valid_priors <- list(loghaz = c("normal","t"),
                         loghr = c("normal","t"),
                         cure = "beta",
                         smooth = "gamma",
                         logor_cure = c("normal","t"),
                         sdnp = "gamma")

    element <- if (is.null(element)) "" else sprintf("[[%s]]",element)
    if (!inherits(prior, "prior")){
        stop(sprintf("`%s%s` should be a call to a prior constructor function such as %s(), see help(survextrap)",
                     priorfullname, element, valid_priors[[priorname]][1]))
    } else {
        if (!(prior$dist %in% valid_priors[[priorname]])){
            vp <- paste0("p_",valid_priors[[priorname]])
            vpstr <- if (length(vp)==1) vp else paste("one of", paste(vp,collapse=","))
            stop(sprintf("`%s%s` should be %s. p_%s was supplied", priorfullname, element, vpstr, prior$dist))
        }
    }
}

get_prior_smooth <- function(prior, est_smooth){
    validate_prior(prior, "smooth")
    if (!est_smooth) return(list(shape=aa(numeric()), rate=aa(numeric())))
    prior
}

## @param prior_list User-supplied prior specification
##
## Can either be like p_normal(0,1) or list(cov1=p_normal(0,2), cov2=p_normal(1, 3))
##
## @param x  List of information describing the linear model , as returned by make_x
get_prior_coveffs <- function(prior, x, modelname){
    priorname <- sprintf("prior_%s", deparse(substitute(prior)))
    if (x$ncovs==0)
        return(list(dist=aa(numeric()),distid=aa(numeric()),
                    location=aa(numeric()), scale=aa(numeric()),
                    df=aa(numeric())))
    else if (is.null(prior))
      prior <- p_normal(0, 2.5)
    prior_list <- validate_prior_bycov(prior, x, priorname)
    for (i in seq_len(x$ncovs)){
      validate_prior(prior_list[[i]], modelname, i)
    }
    dist <- sapply(prior_list, function(x)x$dist)
    distid <- sapply(prior_list, function(x)x$distid)
    location <- sapply(prior_list, function(x)x$location)
    scale <- sapply(prior_list, function(x)x$scale)
    df <- sapply(prior_list, function(x){if (is.null(x$df)) 1 else x$df})
    res <- data.frame(term=x$xnames, dist=dist, distid=distid,
                      location=location, scale=scale, df=df)
    rownames(res) <- NULL
    attr(res, "r") <- lapply(prior_list,function(x)x[c("r")])
    res
}

## Checks on priors specified as a list with one component per covariate
## If supplied as a prior for a single component, replicate that into a list

validate_prior_bycov <- function(prior, x, priorname){
  if (inherits(prior, "prior")){
    prior_list <- rep(list(prior), x$ncovs)
    names(prior_list) <- x$xnames
  } else {
    if (!is.list(prior))
      stop(sprintf("%s should be a call (or list of calls) to a prior constructor function",
                   priorname))
    if (length(prior)!=x$ncovs){
      plural <- if (length(prior)>1) "s" else ""
      stop(sprintf("%s has %s component%s, but there are %s coefficients in the model",
                   priorname, length(prior), plural, x$ncovs))
    }
    if (!identical(sort(names(prior)), sort(x$xnames))){
      quoted_names <- sprintf("\"%s\"", x$xnames)
      stop(sprintf("names of %s do not match names of covariate coefficients: %s",
                   priorname, paste(quoted_names,collapse=",")))
    }
    prior_list <- prior
  }
  prior_list
}

get_prior_sdnp <- function(prior, x, nonprop){
  if (!nonprop) return(array(dim=c(0,2)))
  if (is.null(prior)) prior <- p_gamma(2, 1)
  prior_list <- validate_prior_bycov(prior, x, priorname="prior_sdnp")
  for (i in seq_len(x$ncovs)){
    validate_prior(prior_list[[i]], "sdnp", i)
  }
  shapes <- sapply(prior_list, function(x)x$shape)
  rates <- sapply(prior_list, function(x)x$rate)
  cbind(shapes, rates)
}

prior_basehaz_sample <- function(mspline,
                                 coefs_mean = NULL,
                                 prior_smooth = p_gamma(2,1), # no constants, as they are not allowed in the model
                                 prior_loghaz,
                                 nsim = 100){
  ## process using common code to survextrap
  loghaz_prior <- get_prior_loghaz(prior_loghaz)
  smooth_prior <- get_prior_smooth(prior_smooth, est_smooth=TRUE)
  loghaz <- loghaz_prior$r(nsim)
  smooth <- smooth_prior$r(nsim)

  coefs_mean <- default_coefs_mean(mspline, coefs_mean)
  lcoefs_mean <- log(coefs_mean[-1] / coefs_mean[1])
  np <- length(lcoefs_mean) + 1
  beta <- matrix(nrow=nsim, ncol=np)
  beta[,1] <- 0
  for (j in 2:np){
    beta[,j] <- rlogis(nsim, lcoefs_mean[j-1], smooth)
  }
  coefs <- exp(beta) / rowSums(exp(beta))
  nlist(loghaz, coefs, beta, smooth, haz=exp(loghaz))
}

##' Draw a sample from the prior distribution of the parameters
##' governing the hazard in a survextrap model for given covariates.
##'
##' This is used, for example, in \code{\link{mspline_priorpred}} to visualise the
##' prior distribution around hazard curves implied by a particular M-spline model
##' and parameter priors.
##'
##' Cure and relative survival models are not currently handled.
##'
##' @aliases prior_sample
##'
##' @inheritParams mspline_priorpred
##' @inheritParams survextrap
##'
##' @return A list with components:
##'
##' `loghaz`: Baseline log scale parameter (`log(eta)` in the notation of the manual). For models with covariates, this is at the covariate values supplied in `X`, or at zero if `X` is not supplied.
##'
##' `haz`: Baseline scale parameter (`eta`).
##'
##' `coefs`: Spline coefficients. For non-proportional hazards model with covariates, these are returned at the suppled value of `X`, or at values of zero if `X` is not supplied.
##'
##' `beta`: Multinomial logit-transformed spline coefficients.
##'
##' `smooth`: Smoothing standard deviation for spline coefficients.
##'
##' @name prior_sample
prior_sample <- function(mspline,
                         coefs_mean = NULL,
                         prior_smooth = p_gamma(2,1),
                                        # no constants for the moment as they are not allowed in the model
                         prior_loghaz,
                         prior_loghr = NULL,
                         x = NULL,
                         X = NULL,
                         nonprop = FALSE,
                         prior_sdnp = p_gamma(2,1),
                         nsim = 100){
  sam <- prior_basehaz_sample(mspline,
                              coefs_mean = coefs_mean, prior_smooth = prior_smooth,
                              prior_loghaz = prior_loghaz, nsim=nsim)
  if (is.null(x)) x <- list(ncovs=0, xnames=NULL)
  validate_prior_x(x, X)

  ## baseline HR
  if (x$ncovs > 0){
    loghr_prior <- get_prior_coveffs(prior_loghr, x, "loghr")
    loghr <- lapply(attr(loghr_prior, "r"), function(x)x$r(nsim))
    loghr <- do.call(cbind, loghr)
    X <- unlist(X)
    linpred <- loghr %*% X
    sam$loghaz <- sam$loghaz + linpred
    sam$haz <- exp(sam$loghaz)
  }

  ## nonprop SD
  sdnp_prior <- get_prior_sdnp(prior_sdnp, x, nonprop)
  if (nonprop && !is.null(x)) {
    sdnp <- vector(x$ncovs, mode="list")
    nvars <- ncol(sam$coefs) # number of spline basis terms
    np_linpred <- matrix(0, nrow=nsim, ncol=nvars-1)
    for (i in seq_len(x$ncovs)){
      sdnp[[i]] <- rgamma(nsim, shape=sdnp_prior[i,1], rate=sdnp_prior[i,2])
      b_np <- matrix(rnorm(nsim *(nvars - 1), 0, sdnp[[i]]), nrow=nsim, ncol=nvars-1)
      np_linpred <- np_linpred + b_np * X[i]
    }
    np_linpred <- cbind(0, np_linpred)
    sam$beta <- sam$beta + np_linpred # Can we do this?  Dims?
    sam$coefs <- exp(sam$beta) / rowSums(exp(sam$beta))
  }

  sam
}

##' Summarises the prior for the constant hazard implied by a particular
##' prior on the hazard scale parameter and spline specification.
##'
##' @inheritParams prior_haz
##'
##' @export
prior_haz_const <- function(mspline,
                            prior_loghaz = p_normal(0, 20),
                            nsim = 10000,
                            quantiles = c(0.025, 0.5, 0.975)){
  loghaz_prior <- get_prior_loghaz(prior_loghaz)
  coefs_mean <- mspline_constant_weights(mspline)
  haz_std <-  hsurvmspline(x=mean(mspline$bknots), alpha=1, coefs=coefs_mean,
                           knots=c(mspline$iknots,mspline$bknots), degree=mspline$degree)
  haz <- haz_std * exp(loghaz_prior$q(quantiles))
  names(haz) <- names(quantile(1, quantiles))
  data.frame(haz=haz, mean=rev(1/haz))
}


validate_prior_x <- function(x, X){
  if (x$ncovs==0 && !is.null(X)) {
    warning("Covariate values X supplied, but x$ncovs=0. Ignoring covariates.")
  }
  else {
    if (!is.list(x) ||
        (!("ncovs" %in% names(x))) ||
        (!("xnames" %in% names(x))))
      stop("x should be a list with names `ncovs` and `xnames`")
    if (!is.numeric(x$ncovs)) stop("`x$ncovs` should be numeric")
    if (x$ncovs > 0){
      if (is.null(X))
        stop("x$ncovs > 0, indicating that there are covariates in the model, but covariate values X not supplied")
      if (!is.character(x$xnames)) stop("`x$xnames` should be character")
      if (length(x$xnames) != x$ncovs) stop(sprintf("`x$xnames` is of length %s, but `x$ncovs` is %s",
                                                    length(x$xnames), x$ncovs))
    }
  }
}


##' Compute consequences of priors chosen for a flexible hazard model in
##' survextrap, in terms of a standard deviation measuring the
##' variability over time of the hazard function or the hazard ratio
##' function.
##'
##' The spline model in survextrap allows the hazard to change over time in an
##' arbitrarily flexible manner.  The prior distributions on the parameters of
##' this model have implications for how much we expect the hazard to plausibly
##' vary over time.  These priors are hard to interpret directly, but this
##' function can be used to compute their implications on a more
##' easily-understandable scale.
##'
##' This is done by:
##'
##' (1) simulating a set of parameters from their prior distributions
##'
##' (2) computing the hazard at a fine grid of equally-spaced points spanning
##' the boundary knots
##'
##' (3) calculating the empirical standard deviation of the set of hazards at
##' these points
##'
##' (4) repeatedly performing steps 1-3, and summarising the distribution of the
##' resulting standard deviations.   This is the implied prior for the hazard
##' variability.
##'
##' `prior_haz_sd` computes the SD of the hazard, and the SD of the inverse hazard is also
##' computed.   The inverse hazard at time t is the expected time to the event given survival to t.
##' The hazard ratio between a high and low value (defined by quantiles of values at different times)
##' is also computed. 
##'
##' `prior_hr_sd` computes the SD of the hazard ratio between two covariate values
##' supplied by the user.
##'
##' All of these SDs refer to the variability over time, e.g. a SD of 0 indicates that the
##' hazard (or inverse hazard, or hazard ratio) is constant with time.
##'
##' @inheritParams mspline_args
##' @inheritParams survextrap
##' @inheritParams mspline_priorpred
##'
##' @param hq Quantiles which define the "low" and "high" values to compute the hazard ratio between.
##' By default, this is `c(0.1, 0.9)`, so that the 10\% and 90\% quantiles are used
##' respectively.
##'
##' @param quantiles Quantiles used to summarise the prior predictive distribution
##' of the standard deviation.
##'
##' @name prior_haz
##' @export
prior_haz_sd <- function(mspline,
                         coefs_mean=NULL,
                         prior_smooth=p_gamma(2,1),
                         prior_loghaz = p_normal(0, 20),
                         prior_loghr = NULL,
                         x = NULL,
                         X = NULL,
                         nonprop = FALSE,
                         prior_sdnp = p_gamma(2,1),
                         tmin=0, tmax=10,
                         nsim = 100,
                         hq = c(0.1, 0.9),
                         quantiles = c(0.025, 0.5, 0.975)){
  pred <- mspline_priorpred(iknots=mspline$iknots, bknots=mspline$bknots, degree=mspline$degree,
                               coefs_mean=coefs_mean, prior_smooth=prior_smooth,
                               prior_loghaz=prior_loghaz, prior_loghr=prior_loghr,
                               x=x, X=X, nonprop=nonprop, prior_sdnp=prior_sdnp,
                               tmin=tmin, tmax=tmax,
                            nsim=nsim)
  pred$mean <- 1 / pred$haz
  sd_haz <- tapply(pred$haz, pred$rep, sd)
  sd_mean <- tapply(pred$mean, pred$rep, sd)
  cv_haz <- sd_haz / tapply(pred$haz, pred$rep, mean) # abandoned 
  cv_mean <- sd_mean / tapply(pred$mean, pred$rep, mean)
  hq <- tapply(pred$haz, pred$rep, quantile, probs = hq)
  hq <- do.call("rbind", hq)
  hr <- hq[,2] / hq[,1]
  data.frame(sd_haz = quantile(sd_haz, quantiles),
             sd_mean = quantile(sd_mean, quantiles),
             hr = quantile(hr, quantiles))
}

##' @rdname prior_haz
##' @export
prior_hr_sd <- function(mspline,
                        coefs_mean=NULL,
                        prior_smooth=p_gamma(2,1),
                        prior_loghaz = p_normal(0, 20),
                        prior_loghr = NULL,
                        x = NULL,
                        X = NULL,
                        X0 = NULL,
                        nonprop = FALSE,
                        prior_sdnp = p_gamma(2,1),
                        tmin=0, tmax=10,
                        nsim = 100,
                        quantiles = c(0.025, 0.5, 0.975)){
  pred <- mspline_priorpred(iknots=mspline$iknots, bknots=mspline$bknots, degree=mspline$degree,
                            coefs_mean=coefs_mean, prior_smooth=prior_smooth,
                            prior_loghaz=prior_loghaz, prior_loghr=prior_loghr,
                            x=x, X=X, X0=X0, nonprop=nonprop, prior_sdnp=prior_sdnp,
                            tmin=tmin, tmax=tmax,
                            nsim=nsim)
  sd_hr <- tapply(pred$hr, pred$rep, sd)
  data.frame(sd_hr = quantile(sd_hr, quantiles))
}

default_coefs_mean <- function(mspline, coefs_mean=NULL, logit=FALSE){
  if (is.null(coefs_mean)){
    coefs_mean <- mspline_constant_weights(mspline=mspline, logit=logit)
  } else {
    coefs_mean <- validate_coefs_mean(coefs_mean)
  }
  coefs_mean
}
