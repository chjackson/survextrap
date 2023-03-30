##' Documentation for common M-spline arguments
##'
##' @name mspline_args
##'
##' @param iknots Vector of internal knot locations. If not supplied, \code{df} has to be specified, in which case
#' the default is \code{df - degree - 1} equally spaced knots between the boundary knots.
##' @param bknots Vector with two elements, giving the boundary knot locations
##' @param degree Spline polynomial degree (defaults to 3)
##' @param df Desired number of basis terms, or "degrees of freedom" in the spline.
##' If \code{iknots} is not supplied, the number of internal knots is then chosen to satisfy this.
##'
NULL

mspline_default_iknots <- function(iknots=NULL, bknots, degree, df){
  if (is.null(iknots)) {
    nik <- df - degree  - 1
    iknots <- seq(bknots[1], bknots[2], length.out=nik+2)[-c(1,nik+2)]
  }
  validate_knots(iknots, name="iknots")
  iknots
}

mspline_default <- function(mspline){
  if (!is.list(mspline)) stop("`mspline` should be a list")
  validate_knots(mspline$bknots, name="bknots")
  if (is.null(mspline$degree)) mspline$degree <- 3
  if (is.null(mspline$iknots)){
    if (!is.null(mspline$df)) {
      mspline$iknots <- mspline_default_iknots(bknots=mspline$bknots, degree=mspline$degree, df=mspline$df)
    }
  } else validate_knots(mspline$iknots, name="iknots")
  mspline
}


##' M-spline survival distribution
##'
##' Probability density, distribution, quantile, random generation, hazard,
##' cumulative hazard, mean and restricted mean functions for the
##' M-spline time-to-event model.
##'
##' This can optionally be combined with a cure probability, and / or with
##' a known background hazard trajectory that is a piecewise-constant function of time.
##'
##'
##' @aliases dsurvmspline psurvmspline qsurvmspline rsurvmspline
##' hsurvmspline Hsurvmspline mean_survmspline rmst_survmspline
##'
##' @param x,q,t Vector of times.
##'
##' @param p Vector of probabilities.
##'
##' @param n Number of random numbers to simulate.
##'
##' @param alpha Log scale parameter.
##'
##' @param coefs Spline basis coefficients. These should sum to 1,
##' otherwise they are normalised internally to sum to 1.
##'
##' @param knots Locations of knots on the axis of time, supplied in
##' increasing order.  These include the two boundary knots.
##'
##' In vectorised usage of these functions, the knots and degree must be
##' the same for all alternative times and parameter values.
##'
##' @param degree Spline polynomial degree.
##'
##' @param log,log.p Return log density or probability.
##'
##' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P(X
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
##' Ramsay, J. O. (1988). Monotone regression splines in action. Statistical Science, 3(4): 425–441.
##'
##' Brilleman, S. L., Elci, E. M., Novik, J. B., & Wolfe, R. (2020). Bayesian survival analysis using the rstanarm R package. arXiv preprint arXiv:2002.09633.
##'
##' Wang, W., Yan, J. (2021). Shape-restricted regression splines with R package splines2. Journal of Data Science_, *19*(3), 498-517.
##'
##' @keywords distribution
##'
##' @name Survmspline
NULL

##' Evaluate an M-spline basis matrix at the specified times.
##'
##' Extrapolation beyond the boundary knots is done by assuming that each basis term
##' is constant beyond the boundary at its value at the boundary.
##' This gives a continuous but non-smooth function.   Each basis term is assumed to be
##' zero at times less than zero, since these models are used for hazard functions
##' in survival data.
##'
##' @param times A numeric vector of times at which to evaluate the basis.
##'
##' @param iknots Internal knots
##'
##' @param bknots Boundary knots
##'
##' @param degree Spline degree
##'
##' @param integrate If \code{TRUE}, then the integrated M-spline (I-spline) basis is returned.
##'
##' @return A two-dimensional array.  Rows are the times, and columns are the basis terms.
##'
mspline_basis <- function(times, iknots, bknots, degree=3, integrate = FALSE) {
  validate_knots(iknots, name="iknots")
  validate_knots(bknots, name="bknots")
  tmax <- bknots[2]
  tmin <- bknots[1]
  # evaluate basis at knots first, to set up use of predict()
  basis0 <- splines2::mSpline(iknots, knots = iknots, Boundary.knots = bknots,
                              degree = degree, intercept = TRUE)

  if (integrate) {
    ibasis0 <- splines2::iSpline(iknots, knots = iknots, Boundary.knots = bknots,
                                 degree = degree, intercept = TRUE)
    out <- matrix(nrow=length(times), ncol=ncol(basis0))
    iind <- times <= tmax & times >= tmin
    times_int <- times[iind]
    if (length(times_int) > 0){
        out[iind] <- predict(ibasis0, times_int)
    }
    eind <- which(times > tmax)
    ## Above the upper boundary knot
    times_ext <- times[eind]
    n_ext <- length(times_ext)
    Mmax <- predict(basis0, tmax)
    Imax <- predict(ibasis0, tmax)
    for (i in seq_len(n_ext)){
        out[eind[i],] <- Imax + Mmax*(times_ext[i] - tmax)
    }
    ## Below the lower boundary knot
    eind <- which(times < tmin & times > 0)
    times_ext <- times[eind]
    n_ext <- length(times_ext)
    Mmin <- predict(basis0, tmin)
    Imin <- predict(ibasis0, tmin)
    for (i in seq_len(n_ext)){
        out[eind[i],] <- Mmin*times_int[i]
    }
  } else {
      times <- pmin(times, tmax)
      times <- pmax(times, tmin)
      out <- predict(basis0, times)
  }
  out[times<=0,] <- 0
  aa(out)
}

##' @rdname Survmspline
##' @export
psurvmspline <- function(q, alpha, coefs, knots, degree=3, lower.tail=TRUE, log.p=FALSE,
                         pcure=0, offsetH=0, backhaz=NULL){
  if (!is.null(backhaz)) offsetH <- get_cum_backhaz(q, backhaz)
  if (is.null(pcure)) pcure <- 0
  ind <- att <- NULL
  d <- survmspline_dist_setup(q=q, alpha=alpha, coefs=coefs, knots=knots, pcure=pcure, offsetH=offsetH)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  if (any(ind)){
    log_cumhaz <- Hsurvmspline(x=q, alpha=alpha, coefs=coefs, knots=knots, degree=degree, log=TRUE,
                               offsetH=0)
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
                         pcure=0, offsetH=0, backhaz=NULL){
    if (!is.null(backhaz)) offsetH <- get_cum_backhaz(x, backhaz)
    if (is.null(pcure)) pcure <- 0
    ind <- att <- NULL
    d <- survmspline_dist_setup(q=x, alpha=alpha, coefs=coefs, knots=knots, pcure=pcure, offsetH=offsetH)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    if (any(ind)){
        knots <- sort(knots)
        iknots <- knots[-c(1,length(knots))]
        bknots <- knots[c(1,length(knots))]
        ibasis <- mspline_basis(q, iknots=iknots, bknots=bknots, degree=degree, integrate=TRUE)
        log_cumhaz <- as.vector(alpha) + log(rowSums(coefs * ibasis))
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
                         pcure=0, offseth=0, backhaz=NULL){
    if (!is.null(backhaz)) offseth <- backhaz$hazard[findInterval(x, backhaz$time)]
    if (is.null(pcure)) pcure <- 0
    ind <- att <- NULL
    d <- survmspline_dist_setup(q=x, alpha=alpha, coefs=coefs, knots=knots, pcure=pcure, offseth=offseth)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    if (any(ind)){
        knots <- sort(knots)
        iknots <- knots[-c(1,length(knots))]
        bknots <- knots[c(1,length(knots))]
        basis <- mspline_basis(q, iknots=iknots, bknots=bknots, degree=degree)
        loghaz <- as.vector(alpha) + log(rowSums(coefs * basis))

        pp <- pcure>0
        logdens <- dsurvmspline(x=x[pp], alpha=alpha[pp], coefs=coefs[pp,,drop=FALSE],
                               knots=knots, degree=degree, pcure=0, log=TRUE)
        logsurv <- psurvmspline(q=q[pp], alpha=alpha[pp], coefs=coefs[pp,,drop=FALSE],
                               knots=knots, degree=degree, pcure=pcure[pp], log.p=TRUE)
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
                         pcure=0, offseth=0, offsetH=0, backhaz=NULL){
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
                           log=TRUE, pcure=pcure, offseth=offseth)
    logsurv <- psurvmspline(q, alpha, coefs, knots, degree,
                            log.p=TRUE, lower.tail=FALSE, pcure=pcure, offsetH=offsetH)
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
                         pcure=0, offsetH=0, backhaz=NULL){
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (is.null(pcure)) pcure <- 0
  qgeneric(psurvmspline, p=p, matargs=c("coefs"), scalarargs=c("knots","degree","backhaz"),
           alpha=alpha, coefs=coefs, knots=knots, degree=degree,
           pcure=pcure, offsetH=offsetH, backhaz=backhaz)
}

##' @rdname Survmspline
##' @export
rsurvmspline <- function(n, alpha, coefs, knots, degree=3,
                         pcure=0, offsetH=0, backhaz=NULL){
  if (length(n) > 1) n <- length(n)
  qsurvmspline(p=runif(n), alpha=alpha, coefs=coefs, knots=knots, degree=degree,
               pcure=pcure, offsetH=offsetH, backhaz=backhaz)
}

survmspline_dist_setup <- function(q, alpha, coefs, knots, pcure=0, offsetH=0, offseth=0){
    validate_knots(knots)
    if (!is.matrix(coefs)) coefs <- matrix(coefs, nrow=1)
    lg <- nrow(coefs)
    nret <- max(length(q), nrow(coefs), length(alpha), length(pcure), length(offsetH), length(offseth))
    att <- attributes(q)
    q <- rep(q, length=nret)
    alpha <- rep(alpha, length=nret)
    coefs <- matrix(rep(as.numeric(t(coefs)), length.out = ncol(coefs) * nret),
                    ncol = ncol(coefs), byrow = TRUE)
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
    if (any(!is.numeric(knots))) stop(sprintf("some of `%s` are not numeric", name))
    if (any(knots < 0)) stop(sprintf("some of `%s` are < 0", name))
    if (name=="bknots")
        if (length(knots) != 2) stop("`bknots` should be a vector of length 2")
    ## splines2 handles checking whether internal knots are within boundary
    ## splines2 handles checking the degree
}


##' @rdname Survmspline
##' @export
rmst_survmspline = function(t, alpha, coefs, knots, degree=3, pcure=0, backhaz=NULL){
    if (is.null(pcure)) pcure <- 0
    rmst_generic(psurvmspline, t, start=0,
                 matargs = c("coefs"),
                 unvectorised_args = c("knots","degree","backhaz"),
                 alpha=alpha, coefs=coefs, knots=knots, degree=degree, pcure=pcure, backhaz=backhaz)
}

##' @rdname Survmspline
##' @export
mean_survmspline = function(alpha, coefs, knots, degree=3, pcure=0, backhaz=NULL){
    nt <- if (is.matrix(coefs)) nrow(coefs) else 1
    rmst_generic(psurvmspline, rep(Inf,nt), start=0,
                 matargs = c("coefs"),
                 unvectorised_args = c("knots","degree","backhaz"),
                 alpha=alpha, coefs=coefs, knots=knots, degree=degree, pcure=pcure, backhaz=backhaz)
}

#
# copied from flexsurv
#
rmst_generic <- function(pdist, t, start=0, matargs=NULL, unvectorised_args=NULL, ...)
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
          (1 - pd) / start_p
      }
      res <- try(integrate(fn, start[i], t[i]))
      if (!inherits(res, "try-error"))
          ret[i] <- res$value
  }
  ret[t<start] <- 0
  if (any(is.nan(ret))) warning("NaNs produced")
  ret
}

##' Determine M-spline basis weights which give a constant function.
##'
##' This works by transforming the weights of the corresponding B-spline basis,
##' which are equal if the B-spline is a constant function. 
##'
##' @param mspline A list with components `iknots` (vector of internal knots), `bknots` (vector of boundary knots)
##' and `degree` (polynomial degree) defining an M-spline configuration. 
##'
##' @param logit If \code{TRUE} then the multinomial logit transform of the coefficients
##' is returned.  This is a vector of length one less than the number of coefficients,
##' with the rth element defined by \eqn{log(coefs[r+1] / coefs[1])}.
##'
##' @references Ramsay, J. O. (1988). Monotone regression splines in action. Statistical Science, 3(4): 425–441.
##'
##' @export
mspline_constant_weights <- function(mspline, logit=FALSE){
  mspline <- mspline_default(mspline)
  iknots <- mspline$iknots; bknots <- mspline$bknots; degree <- mspline$degree
  knot_seq <- c(rep(bknots[1], degree+1), iknots, rep(bknots[2], degree+1))
  K <- length(iknots) + degree + 1
  p_const <- numeric(K)
  rescale <- (degree+1)*(bknots[2] - bknots[1])
  for(i in 1:K){
    p_const[i] <- (knot_seq[i+degree+1] - knot_seq[i]) / rescale
  }
  p_const <- p_const/sum(p_const)
  if (logit) log(p_const[-1]/p_const[1]) else p_const
}
