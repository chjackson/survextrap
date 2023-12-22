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
##' Note caution is required about how treatment groups (for example)
##' are stored in your data.  If these are coded as numeric (0/1),
##' then if `newdata` is not specified only one output will be shown,
##' which relates to the average value of this numeric variable over
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
                            newdata0=NULL, wane_period=NULL, wane_nt=10,
                            niter=NULL, summ_fns=NULL, sample=FALSE, ...){
  res <- rmst(x, t=Inf, newdata=newdata, niter=niter, summ_fns=summ_fns,
              newdata0=newdata0, wane_period=wane_period, wane_nt=wane_nt,
              sample=sample)
  res$variable <- "mean"
    res$t <- NULL
  res
}

default_newdata <- function(x, newdata=NULL){
  all_covs <- union(names(x$x$mfbase), names(x$xcure$mfbase))
  all_covs_str <- union(all_covs, x$backhaz_strata)
  has_covs <- (x$ncovs>0) || (x$ncurecovs>0) || (!is.null(x$backhaz_strata))

  if (is.null(newdata) && has_covs){
    newdata <- as.data.frame(c(x$x$mfbase, x$xcure$mfbase))[all_covs]
    str_notcovs <- setdiff(x$backhaz_strata, names(x$x$mfbase)) # strata that are not also covariates
    if (length(str_notcovs) > 0){
      default_strata <- x$backhaz[1,str_notcovs,drop=FALSE]
      if (nrow(newdata)==0)
        newdata <- default_strata
      else
        newdata <- cbind(newdata, default_strata[rep(1,nrow(newdata)),,drop=FALSE])
    }
  }  else {
    if (is.list(newdata)) newdata <- as.data.frame(newdata)
    if (has_covs && !(is.data.frame(newdata)))
      stop("`newdata` should be a data frame or list")
    missing_covs <- all_covs_str[!(all_covs_str %in% names(newdata))]
    if (length(missing_covs) > 0){
      plural <- if (length(missing_covs) > 1) "s" else ""
      stop(sprintf("Values of covariate%s `%s` not included in `newdata`",
                   plural, paste(missing_covs, collapse=",")))
    }
  }
  newdata
}

default_newdata0 <- function(newdata, newdata0=NULL){
  if (is.null(newdata0))
    newdata[rep(1,nrow(newdata)),,drop=FALSE]
  else newdata0
}

##' Restricted mean survival time
##'
##' Compute the restricted mean survival time from a model fitted with
##' \code{\link{survextrap}}.  Defined as the integral of the fitted
##' survival curve up to a specified time.  This relies on numerical
##' integration, which is done for every parameter in the MCMC sample,
##' so it may be slow.
##'
##' @inheritParams mean.survextrap
##'
##' @param t Vector of times.  The restricted mean survival time up to each one of these times
##' will be computed.
##'
##' @examples
##' mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
##' rmst(mod, t=3, niter=100)
##' rmst(mod, t=3, summ_fns=list(mean=mean), niter=100)
##'
##' @return A data frame with each row containing posterior summary statistics
##' for a particular covariate value.
##'
##' An attribute \code{"sample"} is also returned, containing a list of matrices
##' containing samples from the posterior distribution of the RMST.  The list has
##' one component for each row of \code{newdata}, defined as a matrix where the
##' number of rows is the number of MCMC samples, and the number of columns is the
##' number of times in \code{t}.
##'
##' @export
rmst <- function(x,t,newdata=NULL,
                  newdata0=NULL, wane_period=NULL, wane_nt=10,
                  niter=NULL,summ_fns=NULL,sample=FALSE){
  newdata <- default_newdata(x, newdata)
  p <- prepare_pars(x=x, newdata=newdata, t=t, niter=niter, newdata0=newdata0, wane_period=wane_period)

  nd_strata <- make_backhaz_strata(x$backhaz_strata, x$backhaz, newdata)$strata
  res_sam <- array(dim=c(p$nt, p$niter, p$nvals))
  for (j in 1:p$nvals){
    ## for rmst_generic, this should be matrix rather than array 
    coefs1 <- matrix(p$coefs[,,j,], ncol=p$nvars)
    backhaz <- x$backhaz[x$backhaz$stratum == nd_strata[j],] # awkward to vectorise this over nvals
    ## note the offset depends on t so we can't supply a fixed offset
    ## instead have to supply backhaz to enable integration of
    ## psurvmspline over t to get rmst
    if (is.null(wane_period))
      res_sam[,,j] <- rmst_survmspline(t, p$alpha[,,j], coefs1,
                                       x$mspline$knots, x$mspline$degree,
                                       pcure = p$cureprob[,,j],
                                       backhaz = backhaz,
                                       bsmooth=x$mspline$bsmooth)
    else {
      coefs0 <- matrix(p$coefs0[,,j,], ncol=p$nvars)
      res_sam[,,j] <- rmst_survmspline_wane(t, alpha1 = p$alpha[,,j], alpha0 = p$alpha0[,,j],
                                            coefs1 = coefs1, coefs0 = coefs0,
                                            knots = x$mspline$knots,
                                            degree = x$mspline$degree,
                                            pcure1 = p$cureprob[,,j],
                                            pcure0 = p$cureprob0[,,j],
                                            backhaz = backhaz,
                                            bsmooth=x$mspline$bsmooth,
                                            wane_period=wane_period, wane_nt=wane_nt
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
##' @inheritParams rmst
##'
##' @export
irmst <- function(x, t, newdata=NULL, newdata0=NULL, wane_period=NULL, wane_nt=10,
                  niter=NULL, summ_fns=NULL){
    newdata <- default_newdata_comparison(x, newdata)
    rmst1 <- rmst(x, t=t, newdata=newdata[1,,drop=FALSE], newdata0=newdata0[1,,drop=FALSE],
                  wane_period = wane_period, wane_nt=wane_nt, niter=niter)
    rmst2 <- rmst(x, t=t, newdata=newdata[2,,drop=FALSE], newdata0=newdata0[2,,drop=FALSE],
                  wane_period = wane_period, wane_nt=wane_nt, niter=niter)
    irmst_sam <- attr(rmst2, "sample")[[1]] - attr(rmst1, "sample")[[1]]
    irmst_sam <- posterior::as_draws(irmst_sam)
    if (is.null(summ_fns))
      summ_fns <- list(median=median, ~quantile(.x, probs=c(0.025, 0.975)))
    res <- do.call(summary, c(list(irmst_sam), summ_fns))
    res <- cbind(t=as.numeric(res$variable), res)
    res$variable <- "irmst"
    attr(res, "sample") <- irmst_sam
    res
}

newdata_to_X <- function(newdata, x, formula=NULL, xlevs=NULL){
  if (is.null(xlevs)) xlevs <- x$x$xlevs
  if (is.null(formula)) formula <- x$formula
  form <- delete.response(terms(formula))
  if (is.null(newdata))
    X <- matrix(1, nrow=1, ncol=1)
  else {
    X <- model.matrix(form, newdata, xlev=xlevs)
  }
  drop_intercept(X)
}

get_alpha_bycovs <- function(x, stanmat, newdata=NULL){
    if (is.null(newdata) && (x$ncovs > 0)) newdata <- x$x$mfbase
    if (NROW(newdata) > 0){
        X <- newdata_to_X(newdata, x, x$formula, x$x$xlevs) # nvals x npars
    } else X <- NULL
    X <- cbind("(Intercept)"=1, X)
    loghr_names <- grep("loghr", colnames(stanmat), value=TRUE)
    beta <- stanmat[, c("alpha", loghr_names), drop=FALSE]
    beta %*% t(X) # niter x nvals
}

get_cureprob_bycovs <- function(x, stanmat, newdata=NULL){
    if (is.null(newdata)) newdata <- x$xcure$mfbase
    pcure <- if (x$cure) stanmat[, "pcure[1]"] else NULL
    if (x$cure && NROW(newdata) > 0){
        X <- newdata_to_X(newdata, x, x$cure_formula, x$xcure$xlevs)
        X <- cbind("(Intercept)"=1, X)
        logor_names <- grep("logor_cure", colnames(stanmat), value=TRUE)
        logit_intercept <- qlogis(pcure)
        beta <- cbind(logit_intercept, stanmat[, c(logor_names)])
        pcure <- plogis(beta %*% t(X))
    }
    pcure
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

## Form a list of posterior draws that includes intercept parameters for a series of
## user-specific covariate values

get_pars <- function(x, newdata=NULL, niter=NULL){
    stanmat <- get_draws(x)
    if (length(stanmat)==0) stop("Stan model does not contain samples")
    if (is.null(niter))
      niter <- nrow(stanmat)
    else  {
      if (niter > nrow(stanmat)){
        warning(sprintf("niter=%s specified but there are only %s posterior samples.  Using only these %s", niter, nrow(stanmat), nrow(stanmat)))
        niter <- nrow(stanmat)
      }
      stanmat <- stanmat[1:niter,,drop=FALSE]
    }
    gamma <- stanmat[, "gamma[1]", drop=FALSE]
    loghr_names <- sprintf("loghr[%s]",seq(x$ncovs))
    loghr <- if (x$ncovs>0) stanmat[, loghr_names, drop=FALSE] else NULL
    coefs <- get_coefs_bycovs(x, stanmat, newdata)

    ## cure probs for zero covariate values
    pcure <- if (x$cure) stanmat[, "pcure[1]"] else NULL
    ## cure probs for user-specified covariate values
    cureprob_user <- get_cureprob_bycovs(x=x, stanmat=stanmat, newdata=newdata)

    ## log intercept at zero covariate values
    alpha    <- stanmat[, "alpha",  drop = FALSE]
    ## log intercept for user-specified covariate values
    alpha_user <- get_alpha_bycovs(x=x, stanmat=stanmat, newdata=newdata)

    res <- nlist(alpha, alpha_user, gamma, loghr,
                 coefs, pcure, cureprob_user)
    attr(res, "niter") <- niter
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
                                 pcure1=p$cureprob, pcure0=p$cureprob0, offsetH=p$offsetH,
                                 wane_period=wane_period, wane_nt=wane_nt)
  else
    res_sam <- psurvmspline(p$times, p$alpha, p$coefs,
                            x$mspline$knots, x$mspline$degree, lower.tail=FALSE,
                            pcure=p$cureprob, offsetH=p$offsetH,
                            bsmooth=x$mspline$bsmooth)
  res_sam <- array(res_sam, dim=c(p$nt, p$niter, p$nvals)) # can this be functioned?
  res <- summarise_output(res_sam, summ_fns, t, newdata, sample=sample)
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
##' @export
hazard <- function(x, newdata=NULL, t=NULL, tmax=NULL,
                    niter=NULL, summ_fns=NULL, sample=FALSE,
                    newdata0=NULL, wane_period=NULL, wane_nt=10) {
  newdata <- default_newdata(x, newdata)
  if (is.null(t)) t <- default_plottimes(x, tmax)
  p <- prepare_pars(x=x, newdata=newdata, t=t, niter=niter, newdata0=newdata0, wane_period=wane_period)
  if (!is.null(wane_period))
    res_sam <- hsurvmspline_wane(p$times, alpha1=p$alpha, alpha0=p$alpha0,
                                 coefs1=p$coefs, coefs0=p$coefs0,
                                 knots=x$mspline$knots, degree=x$mspline$degree, bsmooth=x$mspline$bsmooth,
                                 pcure1=p$cureprob, pcure0=p$cureprob0, offseth=p$offseth,
                                 wane_period=wane_period, wane_nt=wane_nt)
  else
    res_sam <- hsurvmspline(p$times, p$alpha, p$coefs,
                            x$mspline$knots, x$mspline$degree,
                            pcure=p$cureprob, offseth=p$offseth,
                            bsmooth=x$mspline$bsmooth)
  res_sam <- array(res_sam, dim=c(p$nt, p$niter, p$nvals)) # can this be functioned?
  res <- summarise_output(res_sam, summ_fns, t, newdata, sample=sample)
  res
}

##' Estimates of cumulative hazard from a survextrap model
##'
##' Estimates of the cumulative hazard at particular times, from a \code{\link{survextrap}} model
##'
##' @inheritParams print.survextrap
##' @inheritParams survival
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
                                 pcure1=p$cureprob, pcure0=p$cureprob0, offsetH=p$offsetH,
                                 wane_period=wane_period, wane_nt=wane_nt)
  else
    res_sam <- Hsurvmspline(p$times, p$alpha, p$coefs,
                            x$mspline$knots, x$mspline$degree,
                            pcure=p$cureprob, offsetH=p$offsetH,
                            bsmooth=x$mspline$bsmooth)
  res_sam <- array(res_sam, dim=c(p$nt, p$niter, p$nvals)) # can this be functioned?
  res <- summarise_output(res_sam, summ_fns, t, newdata, sample=sample)
  res
}


default_newdata_comparison <- function(x, newdata){
  if (x$ncovs==0 && x$ncurecovs==0)
    stop("No covariates are in the model")
  if (!is.null(newdata))
    if (!is.data.frame(newdata) || (nrow(newdata) != 2))
      stop("`newdata` should be a data frame with two rows")
  nd <- default_newdata(x, newdata)
  if (!is.null(nd) && is.data.frame(nd) && (nrow(nd) == 2))
    return(nd)
  else {
    stop("`newdata` should be supplied explicitly, as a data frame with two rows, unless the only covariate in the model is a factor with two levels")
  }
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
#' @inheritParams hazard
#'
#' @export
hazard_ratio <- function(x, newdata=NULL, t=NULL, tmax=NULL, niter=NULL, summ_fns=NULL, sample=FALSE) {
  newdata <- default_newdata_comparison(x, newdata)
  if (is.null(t)) t <- default_plottimes(x, tmax)
  haz1 <- hazard(x, newdata=newdata[1,,drop=FALSE], t=t, tmax=tmax, niter=niter, sample=TRUE)[,,1]
  haz2 <- hazard(x, newdata=newdata[2,,drop=FALSE], t=t, tmax=tmax, niter=niter, sample=TRUE)[,,1]
  hr_sam <- haz2 / haz1
  summ_fns <- default_summfns(summ_fns)
  if (!sample){
    res <- apply(hr_sam, 1, do_summfns, summ_fns)
    nsumms <- length(do_summfns(1, summ_fns))
    res <- as.data.frame(t(matrix(res, nrow=nsumms, ncol=length(t))))
    names(res) <- attr(summ_fns, "summnames")
    res <- cbind(t = t, res)
  } else res <- hr_sam
  res
}

deconstruct_mspline <- function(x, scale=1, tmax=NULL){
  sm <- summary(x)
  cf <- sm[sm$variable=="coefs", ]$median
  scale <- exp(sm[sm$variable=="alpha", ]$median)
  bh <- x$mspline
  mspline_plotdata(knots=bh$knots, degree=bh$degree,
                   bsmooth=bh$bsmooth, coefs=cf, scale=scale, tmax=tmax)
}

get_coefs_bycovs <- function(x, stanmat, newdata=NULL, X=NULL){
  if (is.null(X)){
    if (is.null(newdata) && (x$nnphcovs > 0)) newdata <- x$x$mfbase
    if (NROW(newdata) > 0){
      X <- newdata_to_X(newdata, x, x$nph_formula, x$xnph$xlevs) # nvals x nnphcovs
    } else X <- NULL
  }
  niter <- nrow(stanmat)
  nvals <- if (is.null(X)) 1 else nrow(X)

  if (x$nnphcovs > 0){

    b_np_names <- grep("b_np", colnames(stanmat), value=TRUE)
    nvars <- x$mspline$nvars
    b_np <- array(stanmat[,b_np_names], dim=c(niter, x$nnphcovs, nvars-1))

    b_mean <- log(x$coefs_mean[-1] / x$coefs_mean[1])
    b_mean_mat <- matrix(rep(b_mean, niter), nrow=niter, byrow=TRUE)
    b_err_names <- grep("b_err", colnames(stanmat), value=TRUE) # ncovs x nvars
    ssd_names <- grep("ssd", colnames(stanmat), value=TRUE)
    b_err <- array(stanmat[,b_err_names], dim=c(niter, nvars-1))
    ssd <- as.numeric(stanmat[,ssd_names])

    coefs <- array(dim=c(niter, nvals, nvars))
    for (i in 1:nvals){
      locoefs <- array(0, dim=c(niter, nvars-1))
      for (j in 1:x$nnphcovs){
        locoefs <- locoefs + b_mean_mat + b_np[,j,] * X[i,j] + b_err*ssd
      }
      ocoefs <- exp(cbind(0, locoefs))
      coefs[,i,] <- ocoefs / rowSums(ocoefs)
    }
  } else {
    ms_coef_names <- sprintf("coefs[%s]",seq(x$nvars))
    coefs      <- stanmat[, ms_coef_names,  drop = FALSE]
    coefs <- array(coefs, dim=c(niter, 1, x$nvars))[,rep(1,nvals),,drop=FALSE]
  }
  coefs
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
##' @export
hrtime <- function(x, newdata=NULL, niter=NULL, summ_fns=NULL, hq=c(0.1, 0.9)){
  newdata <- default_newdata(x, newdata)
  t <- seq(0, max(x$mspline$knots), length.out=100)
  haz <- hazard(x, newdata=newdata, niter=niter, t=t, sample=TRUE)
  if (!is.numeric(hq) || (length(hq) !=2))
    stop("`hq` should be a numeric vector of length 2")
  haz_hilo <- apply(haz, c(2,3), quantile, hq)
  hr <- haz_hilo[2,,] / haz_hilo[1,,]
  hr <- haz_hilo[2,,] / haz_hilo[1,,]
  dim(hr) <- dim(haz_hilo)[c(2,3)]
  summ_fns <- default_summfns(summ_fns)
  res <- t(apply(hr, 2, do_summfns, summ_fns))
  res <- as.data.frame(res)
  names(res) <- attr(summ_fns, "summnames")
  if (!is.null(newdata))
    res <- cbind(newdata, res)
  res
}

do_summfns <- function(y, summ_fns) {
  dat <- do.call(summary, c(list(posterior::as_draws(matrix(y,ncol=1))),
                            summ_fns))
  vec <- dat[,-1,drop=FALSE]
  setNames(as.numeric(vec), names(dat)[-1])
}

default_summfns <- function(summ_fns){
  if (is.null(summ_fns)){
    summ_fns <- list(~quantile(.x, probs=c(0.5, 0.025, 0.975), na.rm=TRUE))
    attr(summ_fns, "summnames") <- c("median","lower","upper")
  } else attr(summ_fns, "summnames") <- names(do_summfns(1, summ_fns))
  summ_fns
}


#' Organise an array of parameter and input values to supply to
#' vectorised lower-level functions to calculate outputs such
#' as RMST, survival, hazard.
#'
#' @inheritParams mean_survspline
#'
#' @return A list with the following components
#'
#' \code{times} (vector of length nt)
#'
#' \code{alpha}, \code{cureprob}, \code{offseth}, \code{offsetH}
#' arrays of dimension \code{nt x niter x nvals}
#'
#' and \code{coefs}, of dimension \code{nt x niter x nvals x nvars}
#'
#' where
#'
#' \code{nt} is number of time points at which to evaluate the output
#'
#' \code{niter} is the number of MCMC samples
#'
#' \code{nvals} is number of covariate values for which to evaluate the output
#'
#' \code{nvars} is the dimension of the M-spline basis.
#'
#' @noRd
prepare_pars <- function(x, newdata, t, niter=NULL,
                         newdata0=NULL, wane_period=NULL){
    pars <- get_pars(x, newdata=newdata, niter=niter)
    nt <- length(t)
    niter <- attr(pars, "niter")
    nvals <- if (is.null(newdata)) 1 else nrow(newdata)
    nvars <- x$mspline$nvars
    times <- array(rep(t, niter*nvals), dim=c(nt, niter, nvals))
    alpha <- array(rep(pars$alpha_user, each=nt), dim=c(nt, niter, nvals))
    coefs <- array(pars$coefs[rep(1:niter, each=nt),,], dim=c(nt, niter, nvals, nvars))
    cureprob <- if (x$cure) array(rep(pars$cureprob_user, each=nt),
                                  dim = c(nt, niter, nvals)) else NULL

    ## Null covariate values to wane towards in waning models
    if (!is.null(wane_period)){
      newdata0 <- default_newdata0(newdata, newdata0)
      pars0 <- get_pars(x, newdata=newdata0, niter=niter)
      alpha0 <- array(rep(pars0$alpha_user, each=nt), dim=c(nt, niter, nvals))
      coefs0 <- array(pars0$coefs[rep(1:niter, each=nt),,], dim=c(nt, niter, nvals, x$mspline$nvars))
      cureprob0 <- if (x$cure) array(rep(pars0$cureprob_user, each=nt),
                                     dim = c(nt, niter, nvals)) else NULL
    } else alpha0 <- coefs0 <- cureprob0 <- NULL

    ## Offsets for background hazard models
    ## Get stratum for each row of newdata.
    nds <- make_backhaz_strata(x$backhaz_strata, x$backhaz, newdata)
    offseth <- offsetH <- array(0, dim=c(nt, 1, nvals))
    if (!is.null(x$backhaz)){
      for (j in 1:nvals){
        backhaz <- x$backhaz[x$backhaz$stratum == nds$dat[j],]
        offseth[,1,j] <- backhaz$hazard[findInterval(t, backhaz$time)]
        offsetH[,1,j] <- get_cum_backhaz(t, backhaz)
      }
    }
    offsetH <- offsetH[,rep(1,niter),]
    offseth <- offseth[,rep(1,niter),]

    nlist(times,
          alpha, coefs, cureprob,
          alpha0, coefs0, cureprob0,
          offseth, offsetH,
          nt, niter, nvals, nvars)
}

#' Calculate posterior summary statistics of a generic survextrap output
#' function (e.g. RMST, survival, hazard)
#'
#' @param res_sam array of dimension \code{nt}, \code{niter}, \code{nvals},
#' giving the output evaluated at different times, MCMC iterations and
#' covariate values respectively
#'
#' @inheritParams rmst_survspline
#'
#' @return a tibble with \code{nt*nvals} rows, and columns for each
#' posterior summary statistic requested in \code{summ_fns} (typically
#' median or mode and 95\% credible limits).
#'
#' @noRd
summarise_output <- function(res_sam, summ_fns, t, newdata, summ_name=NULL,
                             sample=FALSE){
  if (sample) return(res_sam)
  nt <- dim(res_sam)[1]; niter <- dim(res_sam)[2]; nvals <- dim(res_sam)[3]
  summ_fns <- default_summfns(summ_fns)
  res <- apply(res_sam, c(1,3), do_summfns, summ_fns)
  nsumms <- length(do_summfns(1, summ_fns))
  res <- as.data.frame(matrix(res, nrow=nt*nvals, ncol=nsumms, byrow=TRUE))
  names(res) <- attr(summ_fns, "summnames")
  res <- cbind(t = rep(t, nvals), res)
  if (!is.null(newdata))
    res <- cbind(newdata[rep(1:nvals,each=nt),,drop=FALSE], res)
  if (!is.null(summ_name)){ # TODO why is rmst named but not survival,hazard
    res$variable <- summ_name
    res <- res[c("variable", setdiff(names(res), "variable"))]
  }
  attr(res, "nvals") <- nvals
  tibble::tibble(res)
}
