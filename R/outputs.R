##' Mean survival time
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
##' @param niter Number of MCMC iterations to use to compute credible
##' intervals.  Set to a low value to make this function quicker, at the cost of
##' some approximation error (which may not be important for plotting).
##'
##' @param ... Other options (currently unused).
##'
##' @export
mean.survextrap <- function(x, newdata=NULL, niter=NULL, ...){
    res <- rmst(x, t=Inf, newdata=newdata, niter=niter)
    res$variable <- "mean"
    res$t <- NULL
    res
}

default_newdata <- function(x, newdata){
  all_covs <- union(names(x$x$mfbase), names(x$xcure$mfbase))
  has_covs <- (x$ncovs>0) || (x$ncurecovs>0)
  if (is.null(newdata) && has_covs){
    newdata <- as.data.frame(c(x$x$mfbase, x$xcure$mfbase))[all_covs]
  }  else {
    if (has_covs && !(is.data.frame(newdata)))
      stop("`newdata` should be a data frame")
    missing_covs <- all_covs[!(all_covs %in% names(newdata))]
    if (length(missing_covs) > 0){
      plural <- if (length(missing_covs) > 1) "s" else ""
      stop(sprintf("Values of covariate%s `%s` not included in `newdata`",
                   plural, paste(missing_covs, collapse=",")))
    }
  }
  newdata
}

##' Restricted mean survival time
##'
##' @inheritParams mean.survextrap
##'
##' @param t Vector of times.  The restricted mean survival time up to each one of these times
##' will be computed.
##'
##' @export
rmst <- function(x, t, newdata=NULL, niter=NULL){
    variable <- NULL
    newdata <- default_newdata(x, newdata)
    pars <- get_pars(x, newdata=newdata, niter=niter)
    niter <- attr(pars, "niter")
    nt <- length(t)
    nvals <- if (is.null(newdata)) 1 else nrow(newdata)
    res <- vector(nvals, mode="list")

    for (j in 1:nvals){
        resmat <- matrix(nrow=niter*nvals, ncol=nt)
        colnames(resmat) <- t
        for (i in 1:nt){
          resmat[,i] <- rmst_survmspline(t[i], pars$alpha_user[,j], pars$coefs[,j,],
                                         x$basehaz$knots, x$basehaz$degree,
                                         pcure = pars$cureprob_user, backhaz = x$backhaz)
#          if (x$cure){
#            cureprob_mat <- matrix(rep(pars$cureprob_user, nt), nrow=niter*nvals, ncol=nt)
#            tmat <- matrix(t, nrow=niter*nvals, ncol=nt, byrow=TRUE)
#                resmat <- cureprob_mat*tmat + (1 - cureprob_mat)*resmat
#          }
        }
        sample <- posterior::as_draws(resmat)
        res[[j]] <- summary(sample, median, ~quantile(.x, probs=c(0.025, 0.975)))
        names(res[[j]])[names(res[[j]]) == "variable"] <- "t"
        res[[j]]$t <- as.numeric(res[[j]]$t)
        res[[j]]$variable <- "rmst"
        res[[j]] <- res[[j]][,c("variable", setdiff(names(res[[j]]), "variable"))]
    }
    res <- do.call("rbind", res)
    if (!is.null(newdata))
        res <- cbind(newdata[rep(1:nvals,each=nt),,drop=FALSE], res)
    attr(res, "sample") <- sample
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
  X
}

get_alpha_bycovs <- function(x, stanmat, newdata=NULL){
    if (is.null(newdata) && (x$ncovs > 0)) newdata <- x$x$mfbase
    if (NROW(newdata) > 0){
        X <- newdata_to_X(newdata, x, x$formula, x$x$xlevs) # nvals x npars
        X <- drop_intercept(X)
    } else X <- NULL
    X <- cbind("(Intercept)"=1, X)
    loghr_names <- grep("loghr", colnames(stanmat), value=TRUE)
    beta <- stanmat[, c("alpha", loghr_names), drop=FALSE]
    beta %*% t(X) # niter x nvals
}

get_cureprob_bycovs <- function(x, stanmat, newdata=NULL){
    if (is.null(newdata)) newdata <- x$xcure$mfbase
    pcure <- if (x$cure) stanmat[, "pcure[1]"] else NULL
    if (x$ncurecovs > 0){
        X <- newdata_to_X(newdata, x, x$cure_formula, x$xcure$xlevs)
        X <- drop_intercept(X)
        X <- cbind("(Intercept)"=1, X)
        logor_names <- grep("logor_cure", colnames(stanmat), value=TRUE)
        logit_intercept <- qlogis(pcure)
        beta <- cbind(logit_intercept, stanmat[, c(logor_names)])
        pcure <- plogis(beta %*% t(X))
    }
    pcure
}

### TODO eta is inconsistently named in the code, log scale in stan, natural in vignette

## Form a matrix of posterior draws in "posterior" package draws_matrix format
## Include intercepts at covariate values of zero, as well as at baseline
## (mean for continuous variables and reference levels for factors)

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
    if (is.null(niter)) niter <- nrow(stanmat)
    else stanmat <- stanmat[1:niter,,drop=FALSE]
    gamma <- stanmat[, "gamma[1]", drop=FALSE]
    loghr_names <- sprintf("loghr[%s]",seq(x$ncovs))
    loghr <- if (x$ncovs>0) stanmat[, loghr_names, drop=FALSE] else NULL
    coefs <- get_coefs_bycovs(x, stanmat, newdata)

    ## cure probs for zero covariate values
    pcure <- if (x$cure) stanmat[, "pcure[1]"] else NULL
    ## cure probs for user-specified covariate values
    cureprob_user <- get_cureprob_bycovs(x=x, stanmat=stanmat, newdata=newdata)

    ## log intercept at zero covariate values (including centered factor contrasts)
    alpha    <- stanmat[, "alpha",  drop = FALSE]
    ## log intercept for user-specified covariate values
    alpha_user <- get_alpha_bycovs(x=x, stanmat=stanmat, newdata=newdata)

    res <- nlist(alpha, alpha_user, gamma, loghr, coefs, pcure, cureprob_user)
    attr(res, "niter") <- niter
    res
}

##' Estimates of survival from a survextrap model
##'
##' @inheritParams print.survextrap
##' @inheritParams mean.survextrap
##'
##' @param times Vector of times at which to compute the estimates.
##'
##' @param tmax Maximum time at which to compute the estimates.  If
##' \code{times} is supplied, then this is ignored.  If \code{times} is not
##' supplied, then \code{times} is set to a set of 100 equally spaced
##' time points from 0 to \code{tmax}.  If both \code{tmax} and \code{times}
##' are not supplied, then \code{tmax} is set to the maximum follow up time
##' in the data.
##'
##' @param sample If \code{TRUE} then the MCMC samples are returned instead
##' of being summarised as a median and 95% credible intervals.
##'
##' @export
survival <- function(x, newdata=NULL, times=NULL, tmax=NULL, niter=NULL, sample=FALSE) {
    timedep_output(x, output="survival", newdata=newdata, times=times, tmax=tmax, niter=niter, sample=sample)
}

##' Estimates of hazard from a survextrap model
##'
##' @inheritParams print.survextrap
##' @inheritParams survival
##'
##' @export
hazard <- function(x, newdata=NULL, times=NULL, tmax=NULL, niter=NULL, sample=FALSE) {
  timedep_output(x, output="hazard", newdata=newdata, times=times, tmax=tmax, niter=niter, sample=sample)
}

default_newdata_hr <- function(x, newdata){
  if (x$ncovs==0 && x$ncurecovs==0)
    stop("No covariates are in the model, so can't compute a hazard ratio")
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
#' Intended for use with non-proportional hazards models (\code{survextrap(...,nonprop=TRUE)}).
#'
#' @param newdata A data frame with two rows.  The hazard ratio will be defined as hazard(first row)
#' divided by hazard(second row).
#'
#' @inheritParams hazard
#'
#' @export
hazard_ratio <- function(x, newdata=NULL, times=NULL, tmax=NULL, niter=NULL, sample=FALSE) {
  newdata <- default_newdata_hr(x, newdata)
  if (is.null(times)) times <- default_plottimes(x, tmax)
  haz1 <- hazard(x, newdata=newdata[1,,drop=FALSE], times=times, tmax=tmax, niter=niter, sample=TRUE)[,,1]
  haz2 <- hazard(x, newdata=newdata[2,,drop=FALSE], times=times, tmax=tmax, niter=niter, sample=TRUE)[,,1]
  hr_sam <- haz1 / haz2
  if (!sample){
    res <- apply(hr_sam, 1, function(y)quantile(y, probs=c(0.5, 0.025, 0.975), na.rm=TRUE))
    res <- as.data.frame(t(res))
    names(res) <- c("median","lower","upper")
    res <- cbind(times = times, res)
  } else res <- hr_sam
  res
}

timedep_output <- function(x, output="survival", newdata=NULL, times=NULL, tmax=NULL, niter=NULL, sample=FALSE){
    newdata <- default_newdata(x, newdata)
    pars <- get_pars(x, newdata=newdata, niter=niter)
    if (is.null(times)) times <- default_plottimes(x, tmax)
    nt <- length(times)
    niter <- attr(pars, "niter")
    nvals <- if (is.null(newdata)) 1 else nrow(newdata)
    times_arr <- array(rep(times, niter*nvals), dim=c(nt, niter, nvals))
    alpha_user_arr <- array(rep(pars$alpha_user, each=nt), dim=c(nt, niter, nvals))
    cureprob_arr <- if (x$cure) array(rep(pars$cureprob_user, each=nt), dim = c(nt, niter, nvals)) else NULL

    ## TODO cumulative hazard
    res_sam <- switch(output,
                      "survival" = survival_core(x, pars, times_arr, alpha_user_arr,
                                                 cureprob_arr, niter, nvals, nt),
                      "hazard" = hazard_core(x, pars, times_arr, alpha_user_arr,
                                             cureprob_arr, niter, nvals, nt)
                      )
    dimnames(res_sam) <- list(time=1:nt, iteration=1:niter, value=1:nvals)

    if (!sample){
        res <- apply(res_sam, c(1,3), function(y)quantile(y, probs=c(0.5, 0.025, 0.975), na.rm=TRUE))
        res <- as.data.frame(matrix(res, nrow=nt*nvals, ncol=3, byrow=TRUE))
        names(res) <- c("median","lower","upper")
        res <- cbind(times = rep(times, nvals), res)
        if (!is.null(newdata))
            res <- cbind(newdata[rep(1:nvals,each=nt),,drop=FALSE], res)
    } else res <- res_sam
    attr(res, "nvals") <- nvals
    res
}

survival_core <- function(x, pars, times_arr, alpha_user_arr, cureprob_arr, niter, nvals, nt){
  coefs_mat <- array(pars$coefs[rep(1:niter, each=nt),,], dim=c(nt*niter*nvals, x$basehaz$nvars))
  surv_sam <- psurvmspline(times_arr, alpha_user_arr, coefs_mat,
                           x$basehaz$knots, x$basehaz$degree, lower.tail=FALSE,
                           pcure=cureprob_arr, backhaz=x$backhaz)
#  if (x$cure)  {
#    surv_sam <- cureprob_arr + (1 - cureprob_arr)*surv_sam
#  }
  surv_sam
}

hazard_core <- function(x, pars, times_arr, alpha_user_arr, cureprob_arr, niter, nvals, nt){
  coefs_mat <- array(pars$coefs[rep(1:niter, each=nt),,], dim=c(nt*niter*nvals, x$basehaz$nvars))
#  if (x$cure)
#    logdens_sam <- dsurvmspline(times_arr, alpha_user_arr, coefs_mat,
#                                x$basehaz$knots, x$basehaz$degree, log=TRUE, backhaz=x$backhaz)
  loghaz_sam <- hsurvmspline(times_arr, alpha_user_arr, coefs_mat,
                             x$basehaz$knots, x$basehaz$degree, log=TRUE,
                             pcure=cureprob_arr, backhaz=x$backhaz)
#  if (x$cure)  {
#    loghaz_sam <- log(1 - cureprob_arr) +  logdens_sam -
#      log(survival_core(x, pars, times_arr, alpha_user_arr, cureprob_arr, niter, nvals, nt))
#  }
  haz_sam <- exp(loghaz_sam)
  haz_sam
}

deconstruct_mspline <- function(x, scale=1, tmax=NULL){
  sm <- summary(x)
  cf <- sm[sm$variable=="coefs", ]$median
  scale <- exp(sm[sm$variable=="alpha", ]$median)
  bh <- x$basehaz
  mspline_plotdata(bh$iknots, bh$bknots, bh$df, bh$degree, cf, scale=10, tmax=tmax)
}

get_coefs_bycovs <- function(x, stanmat, newdata=NULL, X=NULL){
  if (is.null(X)){
    if (is.null(newdata) && (x$ncovs > 0)) newdata <- x$x$mfbase
    if (NROW(newdata) > 0){
      X <- newdata_to_X(newdata, x, x$formula, x$x$xlevs) # nvals x ncovs
      X <- drop_intercept(X)
    } else X <- NULL
  }
  niter <- nrow(stanmat)
  nvals <- if (is.null(X)) 1 else nrow(X)

  if (x$nonprop){
    b_np_names <- grep("b_np", colnames(stanmat), value=TRUE)
    nvars <- x$basehaz$nvars
    b_np <- array(stanmat[,b_np_names], dim=c(niter, x$ncovs, nvars-1))

    b_mean <- log(x$coefs_mean[-1] / x$coefs_mean[1])
    b_mean_mat <- matrix(rep(b_mean, niter), nrow=niter, byrow=TRUE)
    b_err_names <- grep("b_err", colnames(stanmat), value=TRUE) # ncovs x nvars
    ssd_names <- grep("ssd", colnames(stanmat), value=TRUE)
    b_err <- array(stanmat[,b_err_names], dim=c(niter, nvars-1))
    ssd <- as.numeric(stanmat[,ssd_names])

    coefs <- array(dim=c(niter, nvals, nvars))
    for (i in 1:nvals){
      locoefs <- array(0, dim=c(niter, nvars-1))
      for (j in 1:x$ncovs){
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
