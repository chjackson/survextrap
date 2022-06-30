##' M-spline survival distribution
##'
##' Probability density, distribution, quantile, random generation, hazard,
##' cumulative hazard, mean and restricted mean functions for the
##' M-spline time-to-event model.
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
##' @return \code{dsurvmspline} gives the density, \code{psurvmspline} gives the
##' distribution function, \code{hsurvmspline} gives the hazard and
##' \code{Hsurvmspline} gives the cumulative hazard, as described in
##'
##' \code{qsurvmspline} gives the quantile function, which is computed by crude
##' numerical inversion (using \code{\link{qgeneric}}).
##'
##' \code{rsurvmspline} generates random survival times by using
##' \code{qsurvmspline} on a sample of uniform random numbers.  Due to the
##' numerical root-finding involved in \code{qsurvmspline}, it is slow compared
##' to typical random number generation functions.
##'
##' @author Christopher Jackson <chris.jackson@@mrc-bsu.cam.ac.uk>
##'
##' @references M-spline paper, rstanarm paper
##'
##' @keywords distribution
##'
##' @name Survmspline
NULL

##' Evaluate an M-spline basis matrix at the specified times.
##'
##' Extrapolation beyond the boundary knots is done by assuming that each basis term
##' is constant beyond the boundary at its value at the boundary.
##' This gives a continuous but non-smooth function.
##'
##' @param times A numeric vector of times .  TODO DOCUMENT negative
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
    times_ext <- times[eind]
    n_ext <- length(times_ext)
    Mmax <- predict(basis0, tmax)
    Imax <- predict(ibasis0, tmax)
    for (i in seq_len(n_ext)){
        out[eind[i],] <- Imax + Mmax*(times_ext[i] - tmax)
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
psurvmspline <- function(q, alpha, coefs, knots, degree=3, lower.tail=TRUE, log.p=FALSE){
    ind <- att <- NULL
    d <- survmspline_dist_setup(q=q, alpha=alpha, coefs=coefs)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    ## TODO validate knots, degree. force to be the same for all q, pars
    if (any(ind)){
        log_cumhaz <- Hsurvmspline(x=q, alpha=alpha, coefs=coefs, knots=knots, degree=degree, log=TRUE)
        log_surv <- as.numeric(-exp(log_cumhaz))  # ie -exp(alpha)*lp
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
Hsurvmspline <- function(x, alpha, coefs, knots, degree=3, log=FALSE){
    ind <- att <- NULL
    d <- survmspline_dist_setup(q=x, alpha=alpha, coefs=coefs)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    if (any(ind)){
        iknots <- knots[-c(1,length(knots))] # TODO consistent interface for this
        bknots <- knots[c(1,length(knots))]
        ibasis <- mspline_basis(q, iknots=iknots, bknots=bknots, degree=degree, integrate=TRUE)
        log_cumhaz <- as.vector(alpha) + log(rowSums(coefs * ibasis))
        if (log)
            ret[ind] <- log_cumhaz
        else
            ret[ind] <- exp(log_cumhaz)
    }
    attributes(ret) <- att
    ret
}

##' @rdname Survmspline
##' @export
hsurvmspline <- function(x, alpha, coefs, knots, degree=3, log=FALSE){
    ind <- att <- NULL
    d <- survmspline_dist_setup(q=x, alpha=alpha, coefs=coefs)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    if (any(ind)){
        iknots <- knots[-c(1,length(knots))]
        bknots <- knots[c(1,length(knots))]
        basis <- mspline_basis(q, iknots=iknots, bknots=bknots, degree=degree)
        loghaz <- as.vector(alpha) + log(rowSums(coefs * basis))
        if (log)
            ret[ind] <- loghaz
        else
            ret[ind] <- exp(loghaz)
    }
    attributes(ret) <- att
    ret
}

##' @rdname Survmspline
##' @export
dsurvmspline <- function(x, alpha, coefs, knots, degree=3, log=FALSE){
    ind <- att <- NULL
    d <- survmspline_dist_setup(q=x, alpha=alpha, coefs=coefs)
    for (i in seq_along(d)) assign(names(d)[i], d[[i]])
    if (any(ind)){
        loghaz <- hsurvmspline(q, alpha, coefs, knots, degree, log=TRUE)
        logsurv <- psurvmspline(q, alpha, coefs, knots, degree, log.p=TRUE, lower.tail=FALSE)
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
qsurvmspline <- function(p, alpha, coefs, knots, degree=3, lower.tail=TRUE, log.p=FALSE){
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    qgeneric(psurvmspline, p=p, matargs=c("coefs"), scalarargs=c("knots","degree"),
             alpha=alpha, coefs=coefs, knots=knots, degree=degree)
}

##' @rdname Survmspline
##' @export
rsurvmspline <- function(n, alpha, coefs, knots, degree=3){
    if (length(n) > 1) n <- length(n)
    ret <- qsurvmspline(p=runif(n), alpha=alpha, coefs=coefs, knots=knots, degree=degree)
    ret
}

survmspline_dist_setup <- function(q, alpha, coefs){
    if (!is.matrix(coefs)) coefs <- matrix(coefs, nrow=1)
    lg <- nrow(coefs)
    nret <- max(length(q), nrow(coefs), length(alpha))
    att <- attributes(q)
    q <- rep(q, length=nret)
    alpha <- rep(alpha, length=nret)
    coefs <- matrix(rep(as.numeric(t(coefs)), length.out = ncol(coefs) * nret),
                    ncol = ncol(coefs), byrow = TRUE)
    ret <- numeric(nret)
    nas <- is.na(q) | is.na(alpha) | is.na(rowSums(coefs))
    ret[nas] <- NA
    nans <- is.nan(q) | is.nan(alpha) | is.nan(rowSums(coefs))
    ret[nans] <- NaN
    ind <- !(nas | nans)
    q <- q[ind]
    alpha <- alpha[ind]
    coefs <- coefs[ind,,drop=FALSE]
    nlist(ret, q, alpha, coefs, ind, att)
}


##' @rdname Survmspline
##' @export
rmst_survmspline = function(t, alpha, coefs, knots, degree=3){
    rmst_generic(psurvmspline, t, start=0,
                 matargs = c("coefs"),
                 scalarargs = c("knots","degree"),
                 alpha=alpha, coefs=coefs, knots=knots, degree=degree)
}

##' @rdname Survmspline
##' @export
mean_survmspline = function(alpha, coefs, knots, degree=3){
    nt <- if (is.matrix(coefs)) nrow(coefs) else 1
    rmst_generic(psurvmspline, rep(Inf,nt), start=0,
                 matargs = c("coefs"),
                 scalarargs = c("knots","degree"),
                 alpha=alpha, coefs=coefs, knots=knots, degree=degree)
}

#
# copied from flexsurv
#
rmst_generic <- function(pdist, t, start=0, matargs=NULL, scalarargs=NULL, ...)
{
  args <- list(...)
  args_mat <- args[matargs]
  args_scalar <- args[scalarargs]
  args[c(matargs,scalarargs)] <- NULL
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
      pdargs <- c(list(start[i]), fargs_vec, fargs_mat, args_scalar)
      start_p <- 1 - do.call(pdist, pdargs)
      fn <- function(end){
          pdargs <- c(list(end), fargs_vec, fargs_mat, args_scalar)
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
