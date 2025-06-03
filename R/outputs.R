##' Mean survival time
##'
##' Compute the mean survival time from a model fitted with
##' \code{\link{survextrap}}.  Defined as the integral of the fitted
##' survival curve from zero to infinity.  This relies on numerical
##' integration, which is done for every parameter in the MCMC sample,
##' so it may be slow.
##'
##' Additionally for some models, the integration up to infinity may
##' not converge, giving an error message.  This typically occurs if
##' there is a substantial probability of high survival times or zero
##' hazards at later times.  The restricted mean survival time can
##' usually be computed in these situations with \code{\link{rmst}},
##' but the model should also be investigated to ensure the posterior
##' distributions are realistic, and simplified or supplemented with
##' external data or informative priors if appropriate.
##'
##' @inheritParams print.survextrap
##'
##' @param newdata Data frame of covariate values to compute the output for.
##' If there are covariates in the model and this is not supplied,
##' the following default is used:
##'
##' (a) if the only covariate is one factor variable, then the output is computed
##' for each level of this factor.
##'
##' (b) if there are multiple covariates, or any numeric covariates, then the output
##' is computed at the mean of each numeric covariate in the original data, and at the
##' baseline level of each factor covariate.
##'
##' Note: caution is required about how treatment groups (for example)
##' are stored in your data.  If these are coded as numeric (0/1),
##' then if `newdata` is not specified, only one output will be shown.
##' This relates to the average value of this numeric variable over
##' the data, which doesn't correspond to either of the treatment
##' groups.  To avoid this, a treatment group should be stored as a
##' factor.
##'
##' @param newdata0 Data frame of covariate values defining the "untreated" group
##' for use in treatment waning models. See \code{\link{Survmspline_wane}}.
##'
##' @param wane_period Vector of two numbers, defining the time period over which
##' the hazard is interpolated between the hazard of the "treated" group (taken from `newdata`)
##' and the hazard of the "untreated" group (taken from `newdata0`).  Optional - if
##' this is not supplied, then no waning is assumed.
##'
##' @param wane_nt Number of intervals defining the piecewise constant approximation
##' to the hazard during the waning period.
##'
##' @param disc_rate Discounting rate used to calculate the discounted mean
##' or restricted mean survival time, using an exponential discounting function.
##'
##' @param niter Number of MCMC iterations to use to compute credible
##' intervals.  Set to a low value to make this function quicker, at the cost of
##' some approximation error (which may not be important for plotting or model
##' development).
##'
##' @param summ_fns A list of functions to use to summarise the posterior sample.
##' This is passed to \code{\link[posterior:summarise_draws]{posterior::summarise_draws}}.
##' By default this is \code{list(median=median, ~quantile(.x, probs=c(0.025, 0.975)))}.
##' If the list is named, then the names will be used for the columns of the
##' output.
##'
##' @param sample If \code{TRUE} then an MCMC sample is returned from the posterior
##' of the output, rather than summary statistics.
##'
##' @param ... Other options (currently unused).
##'
##' @return A data frame with each row containing posterior summary statistics
##' for a particular covariate value.
##'
##' An attribute \code{"sample"} is also returned, containing a matrix
##' of samples from the posterior distribution of the RMST.
##'
##' @export
mean.survextrap <- function(x, newdata=NULL,
                            newdata0=NULL, wane_period=NULL, wane_nt=10, disc_rate = 0,
                            niter=NULL, summ_fns=NULL, sample=FALSE, ...){
  res <- rmst(x, t=Inf, newdata=newdata, niter=niter, summ_fns=summ_fns,
              newdata0=newdata0, wane_period=wane_period, wane_nt=wane_nt,
              sample=sample, disc_rate=disc_rate, method="adaptive")
  res$variable <- "mean"
  res$t <- NULL
  res
}

##' Restricted mean survival time
##'
##' Compute the restricted mean survival time from a model fitted with
##' \code{\link{survextrap}}.  Defined as the integral of the fitted
##' survival curve up to a specified time.
##'
##' @inheritParams mean.survextrap
##'
##' @param t Vector of times.  The restricted mean survival time up to each one of these times
##' will be computed.
##'
##' @param method Method of numerical integration to obtain the
##'   restricted mean survival time from the survival function.
##'
##' The default is `method="gl"`, a Gauss-Legendre method with 100
##' nodes between zero and the maximum of `t`.
##'
##' `method="adaptive"` uses the base R `integrate` function, which is
##' much slower, but potentially more robust for badly-behaved
##' survival functions.
##'
##' @param gl_nodes Number of nodes for the Gauss-Legendre method.
##'
##' @examples
##' mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
##' rmst(mod, t=3, niter=100)
##' rmst(mod, t=3, summ_fns=list(mean=mean), niter=100)
##'
##' @return A data frame (tibble) with each row containing posterior summary statistics
##' for a particular time and covariate value.
##'
##' Or if \code{sample=TRUE}, an array with dimensions
##' \code{length(t)}, \code{niter}, \code{nrow(newdata)}, giving the
##' RMST evaluated at different times, MCMC iterations and covariate
##' values respectively.
##'
##' @export
rmst <- function(x,t,newdata=NULL,
                  newdata0=NULL, wane_period=NULL, wane_nt=10, disc_rate=0, method="gl", gl_nodes=100,
                  niter=NULL,summ_fns=NULL,sample=FALSE){
  if (any(t==Inf)) method <- "adaptive"
  newdata <- default_newdata(x, newdata)

  p <- prepare_pars(x=x, newdata=newdata, niter=niter, newdata0=newdata0, wane_period=wane_period,
                    t = NULL) # result is crossed times x pars inside rmst_generic

  res_sam <- array(dim=c(length(t), p$niter, p$nvals))
  for (j in 1:p$nvals){
    ## for rmst_generic, this should be matrix rather than array
    coefs1 <- matrix(p$coefs[,,j,], ncol=p$nvars)
    if (is.null(wane_period))
      res_sam[,,j] <- rmst_survmspline(t, p$alpha[,,j], coefs1,
                                       x$mspline$knots, x$mspline$degree,
                                       pcure = p$pcure[,,j],
                                       backhaz = p$backhaz[[j]],
                                       bsmooth=x$mspline$bsmooth,
                                       disc_rate = disc_rate,
                                       method = method, gl_nodes=gl_nodes)
    else {
      coefs0 <- matrix(p$coefs0[,,j,], ncol=p$nvars)
      res_sam[,,j] <- rmst_survmspline_wane(t, alpha1 = p$alpha[,,j], alpha0 = p$alpha0[,,j],
                                            coefs1 = coefs1, coefs0 = coefs0,
                                            knots = x$mspline$knots,
                                            degree = x$mspline$degree,
                                            pcure1 = p$pcure[,,j],
                                            pcure0 = p$pcure0[,,j],
                                            backhaz = p$backhaz[[j]],
                                            bsmooth=x$mspline$bsmooth,
                                            wane_period=wane_period, wane_nt=wane_nt,
                                            disc_rate = disc_rate,
                                            method = method, gl_nodes=gl_nodes
                                            )
    }
  }
  res <- summarise_output(res_sam, summ_fns, t, newdata,
                          summ_name="rmst", sample=sample)
  res
}

##' Incremental restricted mean survival time
##'
##' Compute the difference in the restricted mean survival times
##' between two covariate values (e.g. treatment groups).
##'
##' The posterior distribution is obtained by calling
##' \code{\link{rmst}} for each group, obtaining each posterior sample
##' from the \code{"sample"} attribute, and taking the difference to
##' get a posterior sample for the difference.
##'
##' @param newdata A data frame with two rows.  The result will be the
##' restricted mean for the covariates in the second row, minus the
##' restricted mean for the covariates in the first row.  If \code{newdata} is
##' omitted for models where the only covariate is a factor with two
##' levels, then this is taken from these levels.  Otherwise \code{newdata}
##' must be supplied explicitly.
##'
##' @return A data frame (tibble) with each row containing posterior summary statistics
##' for a particular time and covariate value.
##'
##' Or if \code{sample=TRUE}, an array with dimensions
##' \code{length(t)}, \code{niter}, \code{nrow(newdata)}, giving the
##' incremental RMST evaluated at different times, MCMC iterations and covariate
##' values respectively.
##'
##' @inheritParams rmst
##'
##' @export
irmst <- function(x, t, newdata=NULL, newdata0=NULL, wane_period=NULL, wane_nt=10,
                  niter=NULL, summ_fns=NULL, sample=FALSE, disc_rate = 0, method="gl", gl_nodes=100){
    newdata <- default_newdata_comparison(x, newdata)
    rmst1 <- rmst(x, t=t, newdata=newdata[1,,drop=FALSE], newdata0=newdata0[1,,drop=FALSE],
                  wane_period = wane_period, wane_nt=wane_nt, niter=niter, sample=TRUE, disc_rate = disc_rate, method=method, gl_nodes=gl_nodes)
    rmst2 <- rmst(x, t=t, newdata=newdata[2,,drop=FALSE], newdata0=newdata0[2,,drop=FALSE],
                  wane_period = wane_period, wane_nt=wane_nt, niter=niter, sample=TRUE, disc_rate = disc_rate, method=method, gl_nodes=gl_nodes)
    irmst_sam <- rmst2 - rmst1
    res <- summarise_output(irmst_sam, summ_fns, t, newdata=NULL,
                            summ_name="irmst", sample=sample)
    res
}

##' Estimates of survival from a \code{\link{survextrap}} model
##'
##' Estimates of survival probabilities at particular times, from a
##' \code{\link{survextrap}} model
##'
##' @inheritParams print.survextrap
##' @inheritParams mean.survextrap
##'
##' @param t Vector of times at which to compute the estimates.
##'
##' @param tmax Maximum time at which to compute the estimates.  If
##' \code{t} is supplied, then this is ignored.  If \code{t} is not
##' supplied, then \code{t} is set to a set of 100 equally spaced
##' time points from 0 to \code{tmax}.  If both \code{tmax} and \code{t}
##' are not supplied, then \code{tmax} is set to the maximum follow up time
##' in the data.
##'
##' @param sample If \code{TRUE} then the MCMC samples are returned instead
##' of being summarised as a median and 95% credible intervals.
##'
##' @return A data frame (tibble) giving posterior summary statistics,
##'   or (if \code{sample=TRUE}) an array giving samples from the
##'   posterior.
##'
##' @export
survival <- function(x, newdata=NULL, t=NULL, tmax=NULL,
                      niter=NULL, summ_fns=NULL, sample=FALSE,
                      newdata0=NULL, wane_period=NULL, wane_nt=10) {
  newdata <- default_newdata(x, newdata)
  if (is.null(t)) t <- default_plottimes(x, tmax)
  p <- prepare_pars(x=x, newdata=newdata, t=t, niter=niter, newdata0=newdata0, wane_period=wane_period)
  if (!is.null(wane_period))
    res_sam <- psurvmspline_wane(p$times, alpha1=p$alpha, alpha0=p$alpha0,
                                 coefs1=p$coefs, coefs0=p$coefs0,
                                 knots=x$mspline$knots, degree=x$mspline$degree,
                                 bsmooth=x$mspline$bsmooth, lower.tail=FALSE,
                                 pcure1=p$pcure, pcure0=p$pcure0, offsetH=p$offsetH,
                                 wane_period=wane_period, wane_nt=wane_nt)
  else
    res_sam <- psurvmspline(p$times, p$alpha, p$coefs,
                            x$mspline$knots, x$mspline$degree, lower.tail=FALSE,
                            pcure=p$pcure, offsetH=p$offsetH,
                            bsmooth=x$mspline$bsmooth)
  res_sam <- array(res_sam, dim=c(p$nt, p$niter, p$nvals)) # can this be functioned?
  res <- summarise_output(res_sam, summ_fns, t, newdata, sample=sample, summ_name="survival")
  res
}

##' Estimates of hazard from a \code{\link{survextrap}} model
##'
##' Estimates of the hazard function at particular times, from a
##' \code{\link{survextrap}} model
##'
##' @inheritParams print.survextrap
##' @inheritParams survival
##'
##' @return A data frame (tibble) giving posterior summary statistics,
##'   or (if \code{sample=TRUE}) an array giving samples from the
##'   posterior.
##'
##' @export
hazard <- function(x, newdata=NULL, t=NULL, tmax=NULL,
                    niter=NULL, summ_fns=NULL, sample=FALSE,
                    newdata0=NULL, wane_period=NULL, wane_nt=10) {
  newdata <- default_newdata(x, newdata)
  if (is.null(t)) t <- default_plottimes(x, tmax, zero=FALSE)
  p <- prepare_pars(x=x, newdata=newdata, t=t, niter=niter, newdata0=newdata0, wane_period=wane_period)
  if (!is.null(wane_period))
    res_sam <- hsurvmspline_wane(p$times, alpha1=p$alpha, alpha0=p$alpha0,
                                 coefs1=p$coefs, coefs0=p$coefs0,
                                 knots=x$mspline$knots, degree=x$mspline$degree, bsmooth=x$mspline$bsmooth,
                                 pcure1=p$pcure, pcure0=p$pcure0, offseth=p$offseth,
                                 wane_period=wane_period, wane_nt=wane_nt)
  else
    res_sam <- hsurvmspline(p$times, p$alpha, p$coefs,
                            x$mspline$knots, x$mspline$degree,
                            pcure=p$pcure, offseth=p$offseth,
                            bsmooth=x$mspline$bsmooth)
  res_sam <- array(res_sam, dim=c(p$nt, p$niter, p$nvals)) # can this be functioned?
  res <- summarise_output(res_sam, summ_fns, t, newdata, sample=sample, summ_name="hazard")
  res
}

##' Estimates of cumulative hazard from a survextrap model
##'
##' Estimates of the cumulative hazard at particular times, from a \code{\link{survextrap}} model
##'
##' @inheritParams print.survextrap
##' @inheritParams survival
##'
##' @return A data frame (tibble) giving posterior summary statistics,
##'   or (if \code{sample=TRUE}) an array giving samples from the
##'   posterior.
##'
##' @export
cumhaz <- function(x, newdata=NULL, t=NULL, tmax=NULL,
                   niter=NULL, summ_fns=NULL, sample=FALSE,
                   newdata0=NULL, wane_period=NULL, wane_nt=10) {
  newdata <- default_newdata(x, newdata)
  if (is.null(t)) t <- default_plottimes(x, tmax)
  p <- prepare_pars(x=x, newdata=newdata, t=t, niter=niter, newdata0=newdata0, wane_period=wane_period)
  if (!is.null(wane_period))
    res_sam <- Hsurvmspline_wane(p$times, alpha1=p$alpha, alpha0=p$alpha0,
                                 coefs1=p$coefs, coefs0=p$coefs0,
                                 knots=x$mspline$knots, degree=x$mspline$degree, bsmooth=x$mspline$bsmooth,
                                 pcure1=p$pcure, pcure0=p$pcure0, offsetH=p$offsetH,
                                 wane_period=wane_period, wane_nt=wane_nt)
  else
    res_sam <- Hsurvmspline(p$times, p$alpha, p$coefs,
                            x$mspline$knots, x$mspline$degree,
                            pcure=p$pcure, offsetH=p$offsetH,
                            bsmooth=x$mspline$bsmooth)
  res_sam <- array(res_sam, dim=c(p$nt, p$niter, p$nvals)) # can this be functioned?
  res <- summarise_output(res_sam, summ_fns, t, newdata, sample=sample, summ_name="cumhaz")
  res
}



#' Hazard ratio against time in a survextrap model
#'
#' Compute the hazard ratio at a series of time points, estimated from
#' a \code{\link{survextrap}} model.  Intended for use with
#' non-proportional hazards models
#' (\code{survextrap(...,nonprop=TRUE)}).  In proportional hazards
#' models (which \code{\link{survextrap}} fits by default) the hazard
#' ratio is constant with time.
#'
#' @param newdata A data frame with two rows.  The hazard ratio will
#'   be defined as hazard(second row) divided by hazard(first row).
#'   If the only covariate in the model is a factor with two levels,
#'   then \code{newdata} defaults to a data frame containing the
#'   levels of this factor, so that the hazard ratio for the second
#'   level versus the first level is computed.  For any other models,
#'   \code{newdata} must be supplied explicitly.
#'
#' Standardisation (with \code{\link{standardise_to}}) is not supported.
#' This might be done by hand by using \code{hazard(...,sample=TRUE)}
#' to obtain posterior samples for the two standardised hazards separately,
#' then summarising by hand.
#'
#' @inheritParams hazard
#'
#' @return A data frame (tibble) with each row containing posterior summary statistics
#' for different times.
#'
#' Or if \code{sample=TRUE}, an array with dimensions
#' \code{length(t)}, \code{niter}, and 1, giving the
#' incremental RMST evaluated at different times and MCMC iterations
#' respectively.
#'
#' @export
hazard_ratio <- function(x, newdata=NULL, t=NULL, tmax=NULL, niter=NULL, summ_fns=NULL, sample=FALSE) {
  newdata <- default_newdata_comparison(x, newdata)
  if (is.null(t)) t <- default_plottimes(x, tmax, zero=FALSE)
  if (isTRUE(attr(newdata, "std")))
    stop("Standardisation not currently supported in `hazard_ratio`")
  haz1 <- hazard(x, newdata=newdata[1,,drop=FALSE], t=t, tmax=tmax, niter=niter, sample=TRUE)[,,1,drop=FALSE]
  haz2 <- hazard(x, newdata=newdata[2,,drop=FALSE], t=t, tmax=tmax, niter=niter, sample=TRUE)[,,1,drop=FALSE]
  hr_sam <- haz2 / haz1
  summarise_output(hr_sam, summ_fns, t, newdata=NULL, sample=sample, summ_name="hr")
}

##' Hazard ratio between high and low values of the hazard over time
##'
##' This is intended as an intuitive single-number measure of how much
##' a hazard function changes over time.  The hazard is computed on an
##' equally-spaced fine grid between the boundary knots.  The ratio
##' between a "high" and "low" one of these hazard values is computed.
##' For example, if the hazard is constant over time, then this hazard
##' ratio will be 1.
##'
##' @inheritParams hazard
##' @inheritParams prior_haz
##'
##' @return A summary of the posterior distribution of this hazard
##'   ratio from the fitted model, as a data frame with one row per
##'   covariate value requested in `newdata`, and one column for each
##'   posterior summary statistic.
##'
##' Or if \code{sample=TRUE}, an array with dimensions
##' \code{1}, \code{niter}, and \code{nrow(newdata)}, giving the
##' incremental RMST evaluated at different MCMC iterations
##' and covariate values respectively.
##'
##' @export
hrtime <- function(x, newdata=NULL, niter=NULL, summ_fns=NULL, hq=c(0.1, 0.9), sample=FALSE){
  newdata <- default_newdata(x, newdata)
  t <- seq(0, max(x$mspline$knots), length.out=100)
  haz <- hazard(x, newdata=newdata, niter=niter, t=t, sample=TRUE)
  if (!is.numeric(hq) || (length(hq) !=2))
    stop("`hq` should be a numeric vector of length 2")
  haz_hilo <- apply(haz, c(2,3), quantile, hq)
  hr_sam <- haz_hilo[2,,,drop=FALSE] / haz_hilo[1,,,drop=FALSE]
  summarise_output(hr_sam, summ_fns, t=NULL, newdata=newdata, sample=sample, summ_name="hrtime")
}


##' Posterior draws from a survextrap model
##'
##' Return the matrix of draws from the posterior distribution of
##' parameters in a \code{\link{survextrap}} model, with all chains
##' collapsed.
##'
##' @inheritParams mean.survextrap
##'
##' @return A matrix of draws in the \code{draws_matrix} format of the
##'   \pkg{posterior} package.
##'
##' @export
get_draws <- function(x){
  if (x$fit_method %in% c("mcmc","vb"))
    fit <- x$stanfit
  else if (x$fit_method=="opt")
    fit <- x$stanfit$theta_tilde
  else fit <- NULL
  stanmat <- posterior::as_draws_matrix(fit)
  stanmat <- posterior::merge_chains(stanmat)
  stanmat
}

deconstruct_mspline <- function(x, scale=1, tmax=NULL){
  sm <- summary(x)
  cf <- sm[sm$variable=="coefs", ]$median
  scale <- exp(sm[sm$variable=="alpha", ]$median)
  bh <- x$mspline
  mspline_plotdata(knots=bh$knots, degree=bh$degree,
                   bsmooth=bh$bsmooth, coefs=cf, scale=scale, tmax=tmax)
}
