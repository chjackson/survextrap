##' M-spline survival distribution
##'
##' Probability density, distribution, quantile, random generation,
##' hazard, cumulative hazard, mean and restricted mean functions for
##' the M-spline time-to-event model.  This can optionally be combined
##' with a cure probability, and / or with a known background hazard
##' trajectory that is a piecewise-constant function of time.
##'
##'
##' @aliases dsurvmspline psurvmspline qsurvmspline rsurvmspline
##'   hsurvmspline Hsurvmspline mean_survmspline rmst_survmspline
##'
##' @param x,q,t Vector of times.
##'
##' @param p Vector of probabilities.
##'
##' @param n Number of random numbers to simulate.
##'
##' @param alpha Log hazard scale parameter.
##'
##' @param coefs Spline basis coefficients. These should sum to 1,
##'   otherwise they are normalised internally to sum to 1.  Supplied
##'   either as a vector with one element per basis term, or a matrix
##'   with one column per basis term, and rows for alternative values
##'   of the coefficients (in vectorised usage of this function).  If
##'   an array is supplied, it is collapsed into a matrix with number
##'   of columns equal to the final dimension of the array.
##'
##' @param knots Locations of knots on the axis of time, supplied in
##' increasing order.  These include the two boundary knots.
##'
##' In vectorised usage of these functions, the knots and degree must be
##' the same for all alternative times and parameter values.
##'
##' @inheritParams mspline_init
##'
##' @param log,log.p Return log density or probability.
##'
##' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P(X
##' \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}{P(X > x)}.
##'
##' @param pcure Probability of "cure", which defaults to zero.  If this is non-zero, this defines a "mixture cure" version of the M-spline model.
##'
##' @param offseth Constant to be added to the hazard, e.g. representing a "background" risk of death from causes other than the cause of interest.
##'
##' @param offsetH Constant to be added to the cumulative hazard.
##'
##' @param backhaz A data frame that defines the background hazard as a piecewise-constant function of time.
##' This should have two columns, named \code{"time"} and
##'  \code{"hazard"}. Each row gives the "background" hazard between the
##'  specified time and the next time. The first element of
##'  \code{"time"} should be 0, and the final row specifies the hazard at all
##'  times greater than the last element of \code{"time"}.
##'
##' @param disc_rate Discounting rate used to calculate the discounted mean or
##' restricted mean survival time, using an exponential discounting function.
##'
##' @return \code{dsurvmspline} gives the density, \code{psurvmspline} gives the
##' distribution function, \code{hsurvmspline} gives the hazard and
##' \code{Hsurvmspline} gives the cumulative hazard.
##'
##' \code{qsurvmspline} gives the quantile function, which is computed by
##' numerical inversion.
##'
##' \code{rsurvmspline} generates random survival times by using
##' \code{qsurvmspline} on a sample of uniform random numbers.
##'
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##'
##' @details These are the same as the M-splines used to model survival data in `rstanarm`, except that an
##' additional assumption is made that the hazard is constant beyond the boundary knots at its
##' value at the boundary.   This gives a continuous but non-smooth function.
##'
##' The "cure" model can be interpreted in two different ways.  These result in identical probability distribution functions for the event time, hence they are indistinguishable from data:
##'
##' (a) a model where everyone follows the same hazard trajectory that is decreasing through time, with a higher rate of decrease for higher `pcure`.
##'
##' (b) a model where a person either has a constant zero hazard at all times, or a hazard that follows a parametric model (M-spline in this case).  The probability that a person has a zero hazard is `pcure`.
##' This is the "mixture model" interpretation.
##'
##' In the "background hazard" model, the overall hazard is defined as a sum of the background hazard and
##' a cause-specific hazard.   The cause-specific hazard
##'  is modelled with the M-spline model, and the background hazard is assumed
##'  to be a known piecewise-constant function defined by `backhaz`.
##'
##' If both `backhaz` and `pcure` are supplied, then the cure probablity applies only to the cause-specific hazard.
##' That is, the cause-specific hazard decreases to zero, and the overall hazard converges towards
##' the background hazard, as time increases.
##'
##' @references
##'
##' Ramsay, J. O. (1988). Monotone regression splines in action. Statistical Science, 3(4): 425-441.
##'
##' Brilleman, S. L., Elci, E. M., Novik, J. B., & Wolfe, R. (2020). Bayesian survival analysis using the rstanarm R package. arXiv preprint arXiv:2002.09633.
##'
##' Wang, W., Yan, J. (2021). Shape-restricted regression splines with R package splines2. Journal of Data Science_, *19*(3), 498-517.
##'
##' @keywords distribution
##'
##' @name Survmspline
NULL


##' @rdname Survmspline
##' @export
psurvmspline <- function(q, alpha, coefs, knots, degree=3, lower.tail=TRUE, log.p=FALSE,
                         pcure=0, offsetH=0, backhaz=NULL, bsmooth=TRUE){
  if (!is.null(backhaz)) offsetH <- get_cum_backhaz(q, backhaz)
  if (is.null(pcure)) pcure <- 0
  ind <- att <- NULL
  d <- survmspline_dist_setup(q=q, alpha=alpha, coefs=coefs, knots=knots, pcure=pcure, offsetH=offsetH)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  if (any(ind)){
    log_cumhaz <- Hsurvmspline(x=q, alpha=alpha, coefs=coefs, knots=knots, degree=degree, log=TRUE,
                               offsetH=0, bsmooth=bsmooth)
    log_surv <- as.numeric(-exp(log_cumhaz))
    pp <- pcure>0
    log_surv[pp] <- log(pcure[pp] + (1 - pcure[pp])*exp(log_surv[pp]))
    log_surv <- log_surv - offsetH   # log(S) = -H

    log_surv[q==Inf] <- -Inf
    if (log.p && !lower.tail)
      ret[ind] <- log_surv
    else if (!log.p && !lower.tail)
      ret[ind] <- exp(log_surv)
    else if (!log.p && lower.tail)
      ret[ind] <- 1 - exp(log_surv)
    else if (log.p && lower.tail)
      ret[ind] <- log(1 - exp(log_surv))
  }
  attributes(ret) <- att
  ret
}

##' @rdname Survmspline
##' @export
Hsurvmspline <- function(x, alpha, coefs, knots, degree=3, log=FALSE,
                         pcure=0, offsetH=0, backhaz=NULL, bsmooth=TRUE){
  if (!is.null(backhaz)) offsetH <- get_cum_backhaz(x, backhaz)
  if (is.null(pcure)) pcure <- 0
  ind <- att <- NULL
  d <- survmspline_dist_setup(q=x, alpha=alpha, coefs=coefs, knots=knots, pcure=pcure, offsetH=offsetH)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  if (any(ind)){
    knots <- sort(knots)
    ibasis <- mspline_basis(q, knots=knots, degree=degree, integrate=TRUE, bsmooth=bsmooth)
    log_cumhaz <- mspline_sum_basis(ibasis, coefs, as.vector(alpha), log=TRUE)
    pp <- pcure>0
    log_cumhaz[pp] = - (pcure[pp] + (1 - pcure[pp])*(-log_cumhaz[pp])) # since surv = -log_cumhaz
    if (log){
      ret[ind] <- log_cumhaz
      ret[ind][offsetH>0] <- log(exp(ret[ind][offsetH>0]) + offsetH[offsetH>0])
    }
    else {
      ret[ind] <- exp(log_cumhaz)
      ret[ind][offsetH>0] <- ret[ind][offsetH>0] + offsetH[offsetH>0]
    }
  }
  attributes(ret) <- att
  ret
}

##' @rdname Survmspline
##' @export
hsurvmspline <- function(x, alpha, coefs, knots, degree=3, log=FALSE,
                         pcure=0, offseth=0, backhaz=NULL, bsmooth=TRUE){
  if (!is.null(backhaz)) offseth <- backhaz$hazard[findInterval(x, backhaz$time)]
  if (is.null(pcure)) pcure <- 0
  ind <- att <- NULL
  d <- survmspline_dist_setup(q=x, alpha=alpha, coefs=coefs, knots=knots, pcure=pcure, offseth=offseth)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  if (any(ind)){
    knots <- sort(knots)
    basis <- mspline_basis(q, knots=knots, degree=degree, bsmooth=bsmooth)
    loghaz <- mspline_sum_basis(basis, coefs, as.vector(alpha), log=TRUE)
    pp <- pcure>0
    logdens <- dsurvmspline(x=x[pp], alpha=alpha[pp], coefs=coefs[pp,,drop=FALSE],
                            knots=knots, degree=degree, pcure=0, log=TRUE, bsmooth=bsmooth)
    logsurv <- psurvmspline(q=q[pp], alpha=alpha[pp], coefs=coefs[pp,,drop=FALSE],
                            knots=knots, degree=degree, pcure=pcure[pp], log.p=TRUE, lower.tail=FALSE, bsmooth=bsmooth)
    loghaz[pp] <- log(1 - pcure[pp]) + logdens - logsurv

    if (log){
      ret[ind] <- loghaz
      ret[ind][offseth>0] <- log(exp(ret[ind][offseth>0]) + offseth[offseth>0])
    }
    else {
      ret[ind] <- exp(loghaz)
      ret[ind][offseth>0] <- ret[ind][offseth>0] + offseth[offseth>0]
    }
  }
  attributes(ret) <- att
  ret
}

##' @rdname Survmspline
##' @export
dsurvmspline <- function(x, alpha, coefs, knots, degree=3, log=FALSE,
                         pcure=0, offseth=0, offsetH=0, backhaz=NULL,
                         bsmooth=TRUE){
  if (!is.null(backhaz)) {
    offseth <- backhaz$hazard[findInterval(x, backhaz$time)]
    offsetH <- get_cum_backhaz(x, backhaz)
  }
  if (is.null(pcure)) pcure <- 0
  ind <- att <- NULL
  d <- survmspline_dist_setup(q=x, alpha=alpha, coefs=coefs, knots=knots,
                              pcure=pcure, offseth=offseth, offsetH=offsetH)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  if (any(ind)){
    loghaz <- hsurvmspline(q, alpha, coefs, knots, degree,
                           log=TRUE, pcure=pcure, offseth=offseth,
                           bsmooth=bsmooth)
    logsurv <- psurvmspline(q, alpha, coefs, knots, degree,
                            log.p=TRUE, lower.tail=FALSE, pcure=pcure, offsetH=offsetH, bsmooth=bsmooth)
    logdens <- loghaz + logsurv
    if (log)
      ret[ind] <- logdens
    else
      ret[ind] <- exp(logdens)
  }
  attributes(ret) <- att
  ret
}

##' @rdname Survmspline
##' @export
qsurvmspline <- function(p, alpha, coefs, knots, degree=3, lower.tail=TRUE, log.p=FALSE,
                         pcure=0, offsetH=0, backhaz=NULL, bsmooth=TRUE){
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (is.null(pcure)) pcure <- 0
  qgeneric(psurvmspline, p=p, matargs=c("coefs"), scalarargs=c("knots","degree","backhaz","bsmooth"),
           alpha=alpha, coefs=coefs, knots=knots, degree=degree,
           pcure=pcure, offsetH=offsetH, backhaz=backhaz, bsmooth=bsmooth)
}

##' @rdname Survmspline
##' @export
rsurvmspline <- function(n, alpha, coefs, knots, degree=3,
                         pcure=0, offsetH=0, backhaz=NULL, bsmooth=TRUE){
  if (length(n) > 1) n <- length(n)
  qsurvmspline(p=runif(n), alpha=alpha, coefs=coefs, knots=knots, degree=degree,
               pcure=pcure, offsetH=offsetH, backhaz=backhaz, bsmooth=bsmooth)
}

survmspline_dist_setup <- function(q, alpha, coefs, knots, pcure=0, offsetH=0, offseth=0){
  validate_knots(knots)
  if (is.vector(coefs)) coefs <- matrix(coefs, nrow=1)
  if (is.array(coefs) && length(dim(coefs))>2){
    nvars <- dim(coefs)[length(dim(coefs))]
    coefs <- matrix(as.vector(coefs), ncol=nvars)
  }
  lg <- nrow(coefs)
  nret <- max(length(q), nrow(coefs), length(alpha), length(pcure), length(offsetH), length(offseth))
  att <- attributes(q)
  q <- rep(q, length=nret)
  alpha <- rep(alpha, length=nret)
  coefs <- coefs[rep(seq_len(nrow(coefs)), length.out=nret),,drop=FALSE]
  pcure <- rep(pcure, length=nret)
  offsetH <- rep(offsetH, length=nret)
  offseth <- rep(offseth, length=nret)
  ret <- numeric(nret)
  nas <- is.na(q) | is.na(alpha) | is.na(rowSums(coefs)) | is.na(pcure) | is.na(offsetH) | is.na(offseth)
  ret[nas] <- NA
  nans <- is.nan(q) | is.nan(alpha) | is.nan(rowSums(coefs)) | is.nan(pcure) | is.nan(offsetH) | is.nan(offseth)
  ret[nans] <- NaN
  ind <- !(nas | nans)
  q <- q[ind]
  alpha <- alpha[ind]
  coefs <- coefs[ind,,drop=FALSE]
  pcure <- pcure[ind]
  offsetH <- offsetH[ind]
  offseth <- offseth[ind]
  nlist(ret, q, alpha, coefs, ind, pcure, offsetH, offseth, att)
}

validate_knots <- function(knots, name="knots"){
  if (any(!is.numeric(knots)))
    stop(sprintf("All `%s` should be numeric", name))
  if (any(knots <= 0))
    stop(sprintf("All %s should be > 0 ", name))
  ## splines2 handles checking whether internal knots are within boundary
  ## splines2 handles checking the degree
}


##' @rdname Survmspline
##' @export
rmst_survmspline = function(t, alpha, coefs, knots, degree=3, pcure=0, offsetH=0, backhaz=NULL, bsmooth=TRUE, disc_rate = 0){
  if (is.null(pcure)) pcure <- 0
  rmst_generic(psurvmspline, t, start=0,
               matargs = c("coefs"),
               unvectorised_args = c("knots","degree","backhaz","bsmooth"),
               alpha=alpha, coefs=coefs, knots=knots, degree=degree, pcure=pcure, offsetH=offsetH, backhaz=backhaz, bsmooth=bsmooth, disc_rate=disc_rate)
}

##' @rdname Survmspline
##' @export
mean_survmspline = function(alpha, coefs, knots, degree=3, pcure=0, offsetH=0, backhaz=NULL, bsmooth=TRUE, disc_rate=0){
  nt <- if (is.matrix(coefs)) nrow(coefs) else 1
  rmst_generic(psurvmspline, rep(Inf,nt), start=0,
               matargs = c("coefs"),
               unvectorised_args = c("knots","degree","backhaz","bsmooth"),
               alpha=alpha, coefs=coefs, knots=knots, degree=degree, pcure=pcure, offsetH=offsetH, backhaz=backhaz, bsmooth=bsmooth, disc_rate=disc_rate)
}

#
# copied from flexsurv
#
rmst_generic <- function(pdist, t, start=0, matargs=NULL, unvectorised_args=NULL, disc_rate = 0, ...)
{
  args <- list(...)
  args_mat <- args[matargs]
  args_unvectorised <- args[unvectorised_args]
  args[c(matargs,unvectorised_args)] <- NULL
  matlen <- if(is.null(matargs)) NULL else sapply(args_mat, function(x){if(is.matrix(x))nrow(x) else 1})
  veclen <- if (length(args) == 0) NULL else sapply(args, length)
  t_len <- length(t)
  maxlen <- max(c(t_len, veclen, matlen))
  if(length(start) == 1) start <- rep(start, length.out=maxlen)
  na_inds <- rep(FALSE, maxlen)
  for (i in seq(along=args)){
      args[[i]] <- rep(args[[i]], length.out=maxlen)
      na_inds <- na_inds | is.na(args[[i]])
  }
  # exponential discounting function
  if (!is.numeric(disc_rate)) stop("`disc_rate` must be numeric")
  if (length(disc_rate) > 1) stop(sprintf("`disc_rate` must be a scalar. found `length(disc_rate)=`%s",length(disc_rate)))
  disc_fn <- function(t, disc_rate){ exp(-t*disc_rate) }

  t <- rep(t, length.out=maxlen)
  for (i in seq(along=args_mat)){
      if (is.matrix(args_mat[[i]])){
          args_mat[[i]] <- matrix(
              apply(args_mat[[i]], 2, function(x)rep(x, length=maxlen)),
              ncol=ncol(args_mat[[i]]),
              byrow=F
          )
      }
      else args_mat[[i]] <- matrix(args_mat[[i]], nrow=maxlen, ncol=length(args_mat[[i]]), byrow=TRUE)
      na_inds <- na_inds | apply(args_mat[[i]], 1, function(x)any(is.na(x)))
  }
  ret <- numeric(maxlen)
  ret[na_inds] <- NA
  for (i in seq_len(maxlen)[!na_inds]){
      fargs_vec <- lapply(args, function(x)x[i])
      fargs_mat <- lapply(args_mat, function(x)x[i,,drop=FALSE])
      pdargs <- c(list(start[i]), fargs_vec, fargs_mat, args_unvectorised)
      start_p <- 1 - do.call(pdist, pdargs)
      fn <- function(end){
        pdargs <- c(list(end), fargs_vec, fargs_mat, args_unvectorised)
        pd <- do.call(pdist, pdargs)
        (1 - pd)*disc_fn(end, disc_rate = disc_rate) / start_p
      }
      res <- try(integrate(fn, start[i], t[i]))
      if (!inherits(res, "try-error"))
          ret[i] <- res$value
  }
  ret[t<start] <- 0
  if (any(is.nan(ret))) warning("NaNs produced")
  ret
}

## @param basis Matrix giving spline basis
## @param coefs Vector or matrix giving coefficients.  If matrix, it should
## have the same number of rows as coefs
mspline_sum_basis <- function(basis, coefs=NULL, alpha=0, time, log=FALSE) {
  coefs <- validate_coefs(coefs, ncol(basis))
  if (!is.matrix(coefs)) coefs <- rep(coefs, each=nrow(basis))
  if (log)
    haz <- alpha + log(rowSums(basis * coefs))
  else
    haz <- exp(alpha) * rowSums(basis * coefs)
  haz
}

validate_coefs <- function(coefs=NULL, nvars){
  if (is.null(coefs))
    coefs <- rep(1, nvars)
  else {
    nc <- if (is.matrix(coefs)) ncol(coefs) else length(coefs)
    if (nc != nvars)
      stop(sprintf("dimension of `coefs` is %s, should be %s for this M-spline specification", length(coefs), nvars))
  }
  if (is.matrix(coefs))
    coefs <- coefs / rowSums(coefs)
  else
    coefs <- coefs / sum(coefs)
  coefs
}

validate_hscale <- function(hscale){
  if (!is.numeric(hscale)) stop("`hscale` should be numeric")
  if (hscale < 0) stop("`hscale` should be a non-negative number")
}
