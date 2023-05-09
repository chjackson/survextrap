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
#' @return A named list to be used internally by the \code{\link{survextrap}} model fitting functions.
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

#' Derive a normal prior for the log hazard scale parameter based on a guess at survival times
#'
#' Derive a normal prior for the log hazard scale parameter based on a
#' guess at survival times.  The scale parameter is hard to interpret,
#' and depends on the spline knots.  However for any scale parameter,
#' we can determine the spline coefficients that give a constant
#' hazard (\code{\link{mspline_constant_coefs}}).  Therefore if we can
#' guess a typical survival time, we can guess a typical hazard (as 1
#' divided by the survival time) and deduce the scale parameter.  The
#' prior is then constructed by assuming normality on the log scale,
#' and assuming the log upper credible limit is two SDs away from the
#' log median.
#'
#' @inheritParams mspline_constant_coefs
#'
#' @param median Best guess (prior median) for a typical survival time
#'
#' @param upper Upper limit of 95% credible interval for a survival time
#'
#' @return A normal prior in the format returned by \code{\link{p_normal}}, which can
#' be passed directly to the \code{prior_hscale} argument in \code{\link{survextrap}}.
#'
#' @seealso \code{\link{prior_haz_const}}, \code{\link{mspline_constant_coefs}}
#'
#' @export
p_meansurv <- function(median, upper, mspline=NULL){
  coefs_mean <- mspline_constant_coefs(mspline)
  haz_std <-  hsurvmspline(x=mean(mspline$knots), alpha=1, coefs=coefs_mean,
                           knots=mspline$knots, degree=mspline$degree, bsmooth=mspline$bsmooth)
  haz_lower <- 1/upper
  haz_median <- 1/median
  logeta_lower <- log(haz_lower / haz_std)
  logeta_median <- log(haz_median / haz_std)
  logeta_sd <- (logeta_median - logeta_lower) / qnorm(0.975)
  p_normal(logeta_median, logeta_sd)
}

#' Derive a normal prior for the log hazard ratio parameter based on a guess at the hazard ratio
#'
#' Derive a normal prior for the log hazard ratio parameter based on a guess at the hazard ratio.
#' This assumes that the log upper limit is 2 standard deviations away from the log median.
#'
#' @param median Best guess (prior median) for a typical hazard ratio
#'
#' @param upper Upper limit of 95% credible interval for hazard ratio
#'
#' @return A normal prior in the format returned by \code{\link{p_normal}}, which can
#' be passed directly to the \code{prior_loghr} argument in \code{\link{survextrap}}.
#'
#' @export
p_hr <- function(median, upper){
  loghr_median <- log(median)
  loghr_sd <- (log(upper) - log(median)) / qnorm(0.975)
  p_normal(loghr_median, loghr_sd)
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

get_prior_hscale <- function(hscale){
  if (is.null(hscale)) hscale <- p_normal(0, 20) # currently unused: default set in function headers
  else validate_prior(hscale)
  hscale
}

get_priors <- function(hscale, loghr, hsd, cure, logor_cure, x, xcure, est_hsd,
                       nonprop, prior_hrsd){
  hscale <- get_prior_hscale(hscale)
  loghr <- get_prior_coveffs(loghr, x, "loghr")
  hsd <- get_prior_hsd(hsd, est_hsd)
  validate_prior(cure)
  logor_cure <- get_prior_coveffs(logor_cure, xcure, "logor_cure")
  sdnp <- get_prior_hrsd(prior_hrsd, x, nonprop)
  nlist(hscale, loghr, hsd, cure, logor_cure, sdnp)
}

validate_prior <- function(prior, priorname=NULL, element=NULL){
    if (is.null(priorname)) priorname <- deparse(substitute(prior))
    priorfullname <- sprintf("prior_%s", priorname)

    valid_priors <- list(hscale = c("normal","t"),
                         loghr = c("normal","t"),
                         cure = "beta",
                         hsd = "gamma",
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

get_prior_hsd <- function(prior, est_hsd){
    validate_prior(prior, "hsd")
    if (!est_hsd) return(list(shape=aa(numeric()), rate=aa(numeric())))
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

get_prior_hrsd <- function(prior, x, nonprop){
  if (!nonprop) return(array(dim=c(0,2)))
  if (is.null(prior)) prior <- p_gamma(2, 1)
  prior_list <- validate_prior_bycov(prior, x, priorname="prior_hrsd")
  for (i in seq_len(x$ncovs)){
    validate_prior(prior_list[[i]], "sdnp", i)
  }
  shapes <- sapply(prior_list, function(x)x$shape)
  rates <- sapply(prior_list, function(x)x$rate)
  cbind(shapes, rates)
}

prior_basehaz_sample <- function(mspline,
                                 coefs_mean = NULL,
                                 prior_hsd = p_gamma(2,1), # no constants, as they are not allowed in the model
                                 prior_hscale,
                                 nsim = 100){
  ## process using common code to survextrap
  hscale_prior <- get_prior_hscale(prior_hscale)
  hsd_prior <- get_prior_hsd(prior_hsd, est_hsd=TRUE)
  alpha <- hscale_prior$r(nsim)
  hsd <- hsd_prior$r(nsim)

  coefs_mean <- default_coefs_mean(mspline, coefs_mean)
  lcoefs_mean <- log(coefs_mean[-1] / coefs_mean[1])
  np <- length(lcoefs_mean) + 1
  beta <- matrix(nrow=nsim, ncol=np)
  beta[,1] <- 0
  for (j in 2:np){
    beta[,j] <- rlogis(nsim, lcoefs_mean[j-1], hsd)
  }
  coefs <- exp(beta) / rowSums(exp(beta))
  nlist(alpha, coefs, beta, hsd, hscale=exp(alpha))
}

##' Sample from the joint prior of parameters in a survextrap model
##'
##' Draws a sample from the joint prior distribution of the parameters
##' governing a survextrap model for given covariates.  This is used,
##' for example, in \code{\link{prior_sample_hazard}}, to visualise the
##' prior distribution around hazard curves implied by a particular
##' M-spline model and parameter priors.
##'
##' @aliases prior_sample
##'
##' @inheritParams prior_sample_hazard
##' @inheritParams survextrap
##'
##' @return A list with components:
##'
##' `alpha`: Baseline log hazard scale parameter (`log(eta)` in the notation of the manual). For models with covariates, this is at the covariate values supplied in `X`, or at zero if `X` is not supplied.
##'
##' `hscale`: Baseline hazard scale parameter (`eta`).
##'
##' `coefs`: Spline coefficients. For non-proportional hazards model with covariates, these are returned at the suppled value of `X`, or at values of zero if `X` is not supplied.
##'
##' `beta`: Multinomial logit-transformed spline coefficients.
##'
##' `hsd`: Smoothing standard deviation for spline coefficients.
##'
##' If `X0` is supplied, then `alpha0`, `hscale0`, `beta0`, `coefs0` are also returned, representing reference covariate values.
##'
##' `pcure` is returned in cure models (the cure probability).
##'
##' @name prior_sample
prior_sample <- function(mspline,
                         coefs_mean = NULL,
                         prior_hsd = p_gamma(2,1),
                         ## no constants for the moment as they are not allowed in the model
                         prior_hscale,
                         prior_loghr = NULL,
                         X = NULL,
                         X0 = NULL,
                         prior_hrsd = NULL,
                         prior_cure = NULL,
                         nsim = 100){
  sam <- prior_basehaz_sample(mspline,
                              coefs_mean = coefs_mean, prior_hsd = prior_hsd,
                              prior_hscale = prior_hscale, nsim=nsim)

  if (is.null(X)) x <- list(ncovs=0, xnames=NULL)
  else {
    X <- validate_prior_X(X)
    x <- list(ncovs = length(X), xnames = names(X))
  }

  ## baseline hazard scale
  if (x$ncovs > 0){
    loghr_prior <- get_prior_coveffs(prior_loghr, x, "loghr")
    loghr <- lapply(attr(loghr_prior, "r"), function(x)x$r(nsim))
    loghr <- do.call(cbind, loghr[x$xnames])
    if (!is.null(X0)){
      X0 <- validate_prior_X(X0)
      linpred0 <- loghr %*% X0
      sam$alpha0 <- sam$alpha + linpred0
      sam$hscale0 <- exp(sam$alpha0)
    }
    linpred <- loghr %*% X
    sam$alpha <- sam$alpha + linpred
    sam$hscale <- exp(sam$alpha)
  }

  ## nonprop haz models
  if (!is.null(prior_hrsd)) {
    sdnp_prior <- get_prior_hrsd(prior_hrsd, x, nonprop=TRUE)
    sdnp <- vector(x$ncovs, mode="list")
    nvars <- ncol(sam$coefs) # number of spline basis terms
    np_linpred <- matrix(0, nrow=nsim, ncol=nvars-1)
    for (i in seq_len(x$ncovs)){
      sdnp[[i]] <- rgamma(nsim, shape=sdnp_prior[i,1], rate=sdnp_prior[i,2])
      b_np <- matrix(rnorm(nsim *(nvars - 1), 0, sdnp[[i]]),
                     nrow=nsim, ncol=nvars-1)
      if (!is.null(X0))
        np_linpred0 <- np_linpred + b_np * X0[i]
      np_linpred <- np_linpred + b_np * X[i]
    }
    if (!is.null(X0)){
      np_linpred0 <- cbind(0, np_linpred0)
      sam$beta0 <- sam$beta + np_linpred0
      sam$coefs0 <- exp(sam$beta0) / rowSums(exp(sam$beta0))
    }
    np_linpred <- cbind(0, np_linpred)
    sam$beta <- sam$beta + np_linpred
    sam$coefs <- exp(sam$beta) / rowSums(exp(sam$beta))
  }
  if (!is.null(prior_cure)) {
    sam$pcure <- prior_cure$r(nsim)
  }

  sam
}

##' Generate and/or plot a sample from the prior distribution of M-spline hazard curves
##'
##' Generates and/or plots the hazard curves (as functions of time)
##' implied by a prior mean for the spline coefficients (a constant
##' hazard by default) and particular priors for the baseline log
##' hazard and smoothness standard deviation.
##'
##' @aliases prior_sample_hazard plot_prior_hazard
##'
##' @inheritParams survextrap
##' @inheritParams mspline_plotsetup
##'
##' @param nsim Number of simulations to draw
##'
##' @param X Either a list, vector or one-row matrix, with one value
##'   for each covariate in the model.  Defaults to zero.  Factors
##'   must be supplied as multiple numeric (0/1) contrasts.
##'   If there are multiple covariates and a single prior given in
##'   `prior_loghr`, then this prior will be used for all the covariates.
##'   If a list of priors is provided, then both this list and \code{X}
##'   must be named with the names of the covariates. 
##'
##' @param X0 An alternative set of "reference" covariate values
##'   (optional).  The hazard ratio between the hazards at \code{X}
##'   and \code{X0} will also be returned.
##'
##' @return \code{prior_sample_hazard} returns a data frame of the
##'   samples, and \code{plot_prior_hazard} generates a plot.  No
##'   customisation options are provided for the plot function, which
##'   is just intended as a quick check.
##'
##' @name prior_sample_hazard
##' @export
prior_sample_hazard <- function(knots=NULL, df=10, degree=3, bsmooth=TRUE,
                                coefs_mean = NULL,
                                prior_hsd = p_gamma(2,1),
                                prior_hscale = NULL,
                                prior_loghr = NULL,
                                X = NULL,
                                X0 = NULL,
                                prior_hrsd = NULL,
                                tmin=0, tmax=10,
                                nsim=10){
  basis <- mspline_plotsetup(knots=knots, tmin=tmin, tmax=tmax, degree=degree, df=df,
                             bsmooth=bsmooth)
  time <- attr(basis, "times")
  knots <- attr(basis, "knots")
  sam <- prior_sample(mspline = list(knots=knots, degree=degree, bsmooth=bsmooth),
                      coefs_mean=coefs_mean, prior_hsd=prior_hsd,
                      prior_hscale=prior_hscale,
                      prior_hrsd=prior_hrsd,
                      X = X, X0=X0, nsim=nsim)
  hazlist <- vector(nsim, mode="list")
  for (i in 1:nsim){
    haz <- mspline_sum_basis(basis, coefs=sam$coefs[i,],
                             alpha = sam$alpha[i], time=time)
    hazlist[[i]] <- data.frame(haz=haz, time=time)
    hazlist[[i]]$rep <- i
  }
  hazdf <- do.call("rbind", hazlist)
  hazdf$mean <- 1/hazdf$haz
  hazdf$rep <- factor(hazdf$rep)

  if (!is.null(X0)){
    ## comparator hazard using the same simulated parameter values
    ## but different covariate values
    hazlist <- vector(nsim, mode="list")
    for (i in 1:nsim){
      haz <- mspline_sum_basis(basis, coefs=sam$coefs0[i,],
                               alpha = sam$alpha0[i], time=time)
      hazlist[[i]] <- data.frame(haz=haz, time=time)
      hazlist[[i]]$rep <- i
    }
    hazdf0 <- do.call("rbind", hazlist)
    hazdf$haz0 <- hazdf0$haz
    hazdf$hr <- hazdf$haz / hazdf$haz0
  }

  attr(hazdf, "knots") <- knots
  hazdf
}


##' Summarises the prior for the constant hazard implied by a particular
##' prior on the hazard scale parameter and spline specification.
##'
##' Summarises the prior for the constant hazard implied by a particular
##' prior on the hazard scale parameter and M-spline specification,
##' when the spline coefficients are fixed to define a constant hazard
##' using \code{\link{mspline_constant_coefs}}.
##' 
##' @inheritParams prior_haz
##'
##' @seealso \code{\link{p_meansurv}}, \code{\link{mspline_constant_coefs}}
##'
##' @export
prior_haz_const <- function(mspline,
                            prior_hscale = p_normal(0, 20),
                            nsim = 10000,
                            quantiles = c(0.025, 0.5, 0.975)){
  hscale_prior <- get_prior_hscale(prior_hscale)
  coefs_mean <- mspline_constant_coefs(mspline)
  haz_std <-  hsurvmspline(x=mean(mspline$knots), alpha=1, coefs=coefs_mean,
                           knots=mspline$knots, degree=mspline$degree, bsmooth=mspline$bsmooth)
  haz <- haz_std * exp(hscale_prior$q(quantiles))
  names(haz) <- names(quantile(1, quantiles))
  data.frame(haz=haz, mean=rev(1/haz))
}

##' Summarises the prior for the hazard ratio implied by a particular
##' prior on the log hazard ratio
##'
##' Summarises the prior for the hazard ratio implied by a particular
##' prior on the log hazard ratio.   Simply applies an exponential
##' transform to quantiles of the given prior.
##'
##' @inheritParams prior_haz
##'
##' @export
prior_hr <- function(prior_loghr = p_normal(0, 2.5),
                     quantiles = c(0.025, 0.5, 0.975)){
  qh <- qnorm(quantiles, prior_loghr$location, prior_loghr$scale)
  res <- exp(qh)
  names(res) <- format_perc(quantiles)
  res
}

## Copied from internal function in stats package, under GPL.  Copyright (c) R Core.
format_perc <- function (x, digits = max(2L, getOption("digits")), probability = TRUE, 
    use.fC = length(x) < 100, ...) 
{
    if (length(x)) {
        if (probability) 
            x <- 100 * x
        ans <- paste0(if (use.fC) 
            formatC(x, format = "fg", width = 1, digits = digits)
        else format(x, trim = TRUE, digits = digits, ...), "%")
        ans[is.na(x)] <- ""
        ans
    }
    else character(0)
}  

validate_prior_X <- function(X){
  X <- unlist(X)
}


##' Determine priors for time-varying hazards and hazard ratios
##'
##' Computes consequences of priors chosen for the parameters `hsd`
##'  and `hrsd` in a flexible hazard model \code{\link{survextrap}} on
##'  an interpretable scale.  This can be used to calibrate Gamma
##'  priors for these parameters to match interpretable beliefs.
##'
##' The spline model in \code{\link{survextrap}} allows the hazard to
##' change over time in an arbitrarily flexible manner.  The prior
##' distributions on the parameters of this model have implications
##' for how much we expect the hazard to plausibly vary over time.
##' These priors are hard to interpret directly, but this function can
##' be used to compute their implications on a more
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
##' computed.   The inverse hazard at time `t` is the expected time to the event given survival to `t`.
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
##' @inheritParams prior_sample_hazard
##'
##' @param hq Quantiles which define the "low" and "high" values of a
##'   time-varying quantity (hazard in `prior_haz_sd` and the hazard
##'   ratio in `prior_hr_sd`).  The ratio between the high and low
##'   values will be summarised, as a measure of time-dependence.  By
##'   default, this is `c(0.1, 0.9)`, so that the 10% and 90%
##'   quantiles are used respectively.
##'
##' @param quantiles Quantiles used to summarise the implied prior distributions
##' of the simulated quantities.
##'
##' @return A data frame with columns `sd_haz` (SD of the hazard),
##'   `sd_mean` (SD of the inverse hazard) and `hr` (ratio between
##'   high/low hazards) (for \code{\link{prior_haz_sd}}), and rows
##'   giving prior quantiles of these.
##'
##' In \code{\link{prior_hr_sd}}, `sd_hr` is the SD of hazard ratios
##'   over time, and `hrr` is the ratio between high/low hazard ratios.
##' 
##'
##' @name prior_haz
##' @export
prior_haz_sd <- function(mspline,
                         coefs_mean=NULL,
                         prior_hsd=p_gamma(2,1),
                         prior_hscale = p_normal(0, 20),
                         prior_loghr = NULL,
                         X = NULL,
                         prior_hrsd = NULL,
                         tmin=0, tmax=NULL,
                         nsim = 1000,
                         hq = c(0.1, 0.9),
                         quantiles = c(0.025, 0.5, 0.975)){
  pred <- prior_sample_hazard(knots=mspline$knots, degree=mspline$degree,
                              bsmooth=mspline$bsmooth,
                              coefs_mean=coefs_mean, prior_hsd=prior_hsd,
                              prior_hscale=prior_hscale, prior_loghr=prior_loghr,
                              X=X, prior_hrsd=prior_hrsd,
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

## Should this be consistent with hr_survextrap , two rows of X, X0?

##' @rdname prior_haz
##' @export
prior_hr_sd <- function(mspline,
                        coefs_mean=NULL,
                        prior_hsd=p_gamma(2,1),
                        prior_hscale = p_normal(0, 20),
                        prior_loghr = NULL,
                        X,
                        X0,
                        prior_hrsd = NULL,
                        tmin=0, tmax=10,
                        nsim = 100,
                        hq = c(0.1, 0.9),
                        quantiles = c(0.025, 0.5, 0.975)){
  pred <- prior_sample_hazard(knots=mspline$knots, degree=mspline$degree,
                              bsmooth=mspline$bsmooth,
                              coefs_mean=coefs_mean, prior_hsd=prior_hsd,
                              prior_hscale=prior_hscale, prior_loghr=prior_loghr,
                              X=X, X0=X0, prior_hrsd=prior_hrsd,
                              tmin=tmin, tmax=tmax,
                              nsim=nsim)
  sd_hr <- tapply(pred$hr, pred$rep, sd)
  hq <- tapply(pred$hr, pred$rep, quantile, probs = hq)
  hq <- do.call("rbind", hq)
  hrr <- hq[,2] / hq[,1]

  data.frame(sd_hr = quantile(sd_hr, quantiles),
             hrr = quantile(hrr, quantiles))
}

default_coefs_mean <- function(mspline, coefs_mean=NULL, logit=FALSE){
  if (is.null(coefs_mean)){
    coefs_mean <- mspline_constant_coefs(mspline=mspline, logit=logit)
  } else {
    coefs_mean <- validate_coefs_mean(coefs_mean)
  }
  coefs_mean
}
