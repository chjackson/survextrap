### Internal functions used in outputs.R


##' Determine default values for covariates used in survextrap output
##' functions.  Rationale is documented in \code{\link{mean.survextrap}}.
##'
##' @inheritParams mean.survextrap
##'
##' @return Data frame with default covariate values
##'
##' @noRd
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

##' Determine default values for covariates used for the untreated
##' group in waning models.
##'
##' This is taken from the first row of default is the set to the
##' first row of \code{newdata}.  Typically this would be a "null"
##' value, such as the control group in a trial.
##'
##' @inheritParams mean.survextrap
##'
##' @return Data frame of covariate values
##'
##' @noRd
default_newdata0 <- function(newdata, newdata0=NULL){
  if (is.null(newdata0))
    newdata[rep(1,nrow(newdata)),,drop=FALSE]
  else newdata0
}

##' Determine default covariate values to use for comparisons
##'
##' This is used in, e.g. \code{\link{irmst}} and \code{\link{hazard_ratio}}.
##'
##' @inheritParams mean.survextrap
##'
##' @return A data frame of two rows, the second row giving the value
##'   of interest, and the first row giving the comparator value.
##'
##' @noRd
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


#' Organise an array of parameter and input values to supply to
#' vectorised lower-level functions to calculate outputs such
#' as RMST, survival, hazard.
#'
#' @inheritParams mean_survspline
#'
#' @return A list with the following components
#'
#' \code{times} (vector of length \code{nt})
#'
#' \code{alpha}, \code{pcure}, \code{offseth}, \code{offsetH},
#' arrays of dimension \code{nt x niter x nvals}
#'
#' \code{coefs}, array of dimension \code{nt x niter x nvals x nvars}
#'
#' \code{nt}, number of time points at which to evaluate the output
#'
#' \code{niter}, number of MCMC samples
#'
#' \code{nvals}, number of covariate values for which to evaluate the output
#'
#' \code{nvars}, dimension of the M-spline basis.
#'
#' \code{alpha0}, \code{coefs0} and \code{pcure0} are also included
#' for waning models.
#'
#' @noRd
prepare_pars <- function(x, newdata, t, niter=NULL,
                         newdata0=NULL, wane_period=NULL){
    pars <- get_pars_bycovs(x, newdata=newdata, niter=niter)
    nt <- length(t)
    niter <- pars$niter
    nvals <- if (is.null(newdata)) 1 else nrow(newdata)
    nvars <- x$mspline$nvars
    times <- array(rep(t, niter*nvals), dim=c(nt, niter, nvals))
    alpha <- array(rep(pars$alpha, each=nt), dim=c(nt, niter, nvals))
    coefs <- array(pars$coefs[rep(1:niter, each=nt),,], dim=c(nt, niter, nvals, nvars))
    pcure <- if (x$cure) array(rep(pars$pcure, each=nt),
                                  dim = c(nt, niter, nvals)) else NULL

    ## Null covariate values to wane towards in waning models
    if (!is.null(wane_period)){
      newdata0 <- default_newdata0(newdata, newdata0)
      pars0 <- get_pars_bycovs(x, newdata=newdata0, niter=niter)
      alpha0 <- array(rep(pars0$alpha, each=nt), dim=c(nt, niter, nvals))
      coefs0 <- array(pars0$coefs[rep(1:niter, each=nt),,], dim=c(nt, niter, nvals, x$mspline$nvars))
      pcure0 <- if (x$cure) array(rep(pars0$pcure, each=nt),
                                     dim = c(nt, niter, nvals)) else NULL
    } else alpha0 <- coefs0 <- pcure0 <- NULL

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
          alpha, coefs, pcure,
          alpha0, coefs0, pcure0,
          offseth, offsetH,
          nt, niter, nvals, nvars)
}


##' Extract samples of parameters from Stan output and combine with
##' covariate values
##'
##' For each parameter value sampled from the posterior, an equivalent
##' set of parameters is produced for each of the user-supplied
##' covariate values, generally supplied in an argument called
##' \code{newdata}.
##'
##' @inheritParams mean.survmspline
##'
##' @return A list with components:
##'
##' \code{alpha} Log hazard scale at the series of covariate values
##' supplied in \code{newdata}.  Matrix with \code{niter} rows and
##' \code{nvals} columns, where \code{niter} is the number of MCMC
##' iterations requested by the user in the \code{niter} argument to
##' the output function, and \code{nvals} is the number of rows of
##' \code{newdata}, or 1 in models with no covariates.
##'
##' \code{coefs} Spline basis coefficients.  Array with dimensions
##' \code{niter} x \code{nvals} x \code{nvars}, where \code{nvars} is
##' the basis dimension.
##'
##' \code{pcure} Probability of cure in cure models.  Matrix with
##' \code{niter} rows and \code{nvals} columns.
##'
##' \code{niter} Number of iterations requested for calculating
##' outputs.  If this function was called with \code{niter=NULL} then
##' the full number of MCMC iterations performed will be used.  Note
##' the user may not want to use this full number of iterations
##' available for outputs, e.g. for rough plots.
##'
##' @noRd
get_pars_bycovs <- function(x, newdata=NULL, niter=NULL){
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

    ## parameters at user-specified covariate values
    alpha <- get_alpha_bycovs(x=x, stanmat=stanmat, newdata=newdata)
    coefs <- get_coefs_bycovs(x, stanmat, newdata)
    pcure <- get_pcure_bycovs(x=x, stanmat=stanmat, newdata=newdata)

    res <- nlist(alpha, coefs, pcure,
                 niter) # gamma, loghr, _base versions not needed
    res
}

##' Evaluate log hazard scale at user-supplied covariate values
##'
##' @return Matrix with \code{niter} rows and \code{nvals} columns.
##'
##' @noRd
get_alpha_bycovs <- function(x, stanmat, newdata=NULL, X=NULL){
  if (is.null(X)){
    if (is.null(newdata) && (x$ncovs > 0)) newdata <- x$x$mfbase
    if (NROW(newdata) > 0){
        X <- newdata_to_X(newdata, x, x$formula, x$x$xlevs) # nvals x npars
    }
  } else if (nrow(X)==0) return(NULL)
  X <- cbind("(Intercept)"=1, X)
  loghr_names <- grep("loghr", colnames(stanmat), value=TRUE)
  beta <- stanmat[, c("alpha", loghr_names), drop=FALSE]
  beta %*% t(X) # niter x nvals
}

##' Evaluate basis coefficients at user-supplied covariate values in
##' non-proportional hazards models
##'
##' @return Array with dimensions \code{niter} x \code{nvals} x
##'   \code{nvars}, where \code{nvars} is the basis dimension.
##'
##' @noRd
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

##' Evaluate cure probability at user-supplied covariate values
##'
##' @return Matrix with \code{niter} rows and one column for each
##' covariate value.
##'
##' @noRd
get_pcure_bycovs <- function(x, stanmat, newdata=NULL, X=NULL){
    if (is.null(newdata)) newdata <- x$xcure$mfbase
    pcure <- if (x$cure) stanmat[, "pcure[1]"] else NULL
    if (x$cure){
      if (is.null(X) && (NROW(newdata) > 0)){
        X <- newdata_to_X(newdata, x, x$cure_formula, x$xcure$xlevs)
      }
      X <- cbind("(Intercept)"=1, X)
      logor_names <- grep("logor_cure", colnames(stanmat), value=TRUE)
      logit_intercept <- qlogis(pcure)
      beta <- cbind(logit_intercept, stanmat[, c(logor_names)])
      pcure <- plogis(beta %*% t(X))
    }
    pcure
}

##' Convert user-supplied covariate values to a design matrix
##'
##' @noRd
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
  if (!is.null(summ_name)){
    res$variable <- summ_name
    res <- res[c("variable", setdiff(names(res), "variable"))]
  }
  attr(res, "nvals") <- nvals
  tibble::tibble(res)
}

##' Return default functions used to summarise the posterior
##'
##' These are the median and equal-tailed 95% credible limits
##'
##' @return List of functions in the syntax used by the
##'   \code{posterior} package.  An attribute is returned giving names
##'   to use for the output displayed to the user.
##'
##' @noRd
default_summfns <- function(summ_fns){
  if (is.null(summ_fns)){
    summ_fns <- list(~quantile(.x, probs=c(0.5, 0.025, 0.975), na.rm=TRUE))
    attr(summ_fns, "summnames") <- c("median","lower","upper")
  } else attr(summ_fns, "summnames") <- names(do_summfns(1, summ_fns))
  summ_fns
}

##' Evaluate posterior summary functions
##'
##' @param y Vector of MCMC draws
##'
##' @param summ_fns List of functions in the syntax used by the
##'   \code{posterior} package
##'
##' @return A vector with one element for the output of each summary
##'   function, e.g. mean, median.
##'
##' @noRd
do_summfns <- function(y, summ_fns) {
  dat <- do.call(summary, c(list(posterior::as_draws(matrix(y,ncol=1))),
                            summ_fns))
  vec <- dat[,-1,drop=FALSE]
  setNames(as.numeric(vec), names(dat)[-1])
}


