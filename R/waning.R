##' M-spline survival distribution under treatment effect waning
##'
##' This defines the CDF, cumulative hazard and hazard of a survival
##' distribution defined by combining the hazards of two different
##' groups (e.g. "treated" and "untreated") each defined by a standard
##' M-spline model. The log hazards of one group and the other are
##' interpolated over a defined period of time.  This may be used for
##' models where the treatment effect wanes, over a period of time,
##' between an estimated hazard ratio and zero.
##'
##' This distribution is defined as follows.
##'
##' * Between time 0 and `wane_period[1]`, the "treated" hazard is
##' used, as defined by an M-spline with intercept `alpha1`.
##'
##' * Between `wane_period[1]` and `wane_period[2]`, the log hazard is
##' defined by linear interpolation.  The waning period is divided
##' into a number of discrete pieces in which the hazard is assumed to
##' be constant, defined by the hazard at the start of the piece.
##' These hazard values are obtained from the spline model, using an
##' intercept parameter `alpha` (log scale parameter) defined by
##' linearly interpolating between `alpha1` and `alpha0` over the
##' waning period.  The cumulative hazard at any time can then be
##' deduced by adding up contributions on each piece.
##'
##' * After `wane_period[2]`, the "untreated" hazard is used, as defined by an
##' M-spline with intercept `alpha0`.
##'
##' See the [methods vignette](https://chjackson.github.io/survextrap/articles/methods.html) for more details and examples.
##'
##' This can be used to predict the hazard of a person treated with a
##' treatment whose short-term effect is estimated from shorter-term
##' data, but we wish to extrapolate this model over a longer period
##' where the effect is assumed to diminish.
##'
##' This may not work if the hazard is zero or infinite at any point
##' in the waning period (thus the log hazard is indeterminate).  This
##' might typically happen at time 0.
##'
##' @aliases psurvmspline_wane Hsurvmspline_wane hsurvmspline_wane
##' dsurvmspline_wane qsurvmspline_wane rsurvmspline_wane
##' rmst_survmspline_wane mean_survmspline_wane
##'
##' @inheritParams Survmspline
##'
##' @param x,q,t Vector of times.
##'
##' @param alpha1,coefs1,pcure1 log hazard intercept, spline coefficients and cure probability before the start of the waning period ("treated")
##'
##' @param alpha0,coefs0,pcure0 log hazard intercept, spline coefficients and cure probability after the end of the waning period ("untreated")
##'
##' @param wane_period vector of two components giving start and stop of waning
##'   period
##'
##' @param wane_nt time resolution for piecewise constant hazard approximation in the
##'   waning period.  If this is not specified, defaults to dividing the waning period into 10 pieces.
##'
##'
##' @return \code{psurvmspline_wane} gives the CDF, \code{Hsurvmspline_wane} gives the cumulative
##' hazard, \code{hsurvmspline_wane} gives the hazard, \code{dsurvmspline_wane} gives the PDF,
##' \code{qsurvmspline_wane} gives the quantiles, and \code{rsurvmspline_wane} generates random
##' numbers from the distribution.
##'
##' @keywords distribution
##'
##' @name Survmspline_wane
NULL

##' @rdname Survmspline_wane
##' @export
Hsurvmspline_wane <- function(x, alpha1, alpha0, coefs1, coefs0, knots, degree=3,
                              wane_period, wane_nt=10, pcure1=0, pcure0=0,
                              offsetH=0,
                              backhaz=NULL, bsmooth=TRUE, log=FALSE){

  if (!is.null(backhaz)) offsetH <- get_cum_backhaz(x, backhaz)
  if (is.null(pcure1)) pcure1 <- 0
  if (is.null(pcure0)) pcure0 <- 0
  if (length(alpha0) != length(alpha1)) stop("lengths of `alpha0` and `alpha1` should be the same")
  att <- NULL
  d <- survmspline_dist_setup(q=x, alpha=alpha1, coefs=coefs1, knots=knots, pcure=pcure1, offsetH=offsetH)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  x <- q
  alpha1 <- d$alpha; coefs1 <- d$coefs; pcure1 <- d$pcure
  d0 <- survmspline_dist_setup(q=x, alpha=alpha0, coefs=coefs0, knots=knots, pcure=pcure0, offsetH=offsetH)
  alpha0 <- d0$alpha; coefs0 <- d0$coefs; pcure0 <- d0$pcure

  if (wane_period[1] >= wane_period[2]) stop("wane_period[2] should be greater than wane_period[1]")
  dt <- (wane_period[2] - wane_period[1]) / wane_nt

  ret <- numeric(length(x))
  attributes(ret) <- att
  ## times before the waning period
  early <- x <= wane_period[1]
  ret[early] <- Hsurvmspline(x[early], alpha=alpha1[early], coefs=coefs1[early,,drop=FALSE], pcure=pcure1[early],
                             knots=knots, degree=degree, bsmooth=bsmooth)

  if (any(x > wane_period[1])){
    times <- seq(wane_period[1], wane_period[2], by=dt)
    ints <- findInterval(x, times)
    nints <- max(ints)

    ## pre-compute matrix of cumulative hazards H for every cut-point in the waning period
    ## then pick out required one for each x (inefficient perhaps)
    H <- h <- matrix(nrow=length(alpha1), ncol=nints)
    H[,1] <- Hsurvmspline(wane_period[1], alpha1, coefs=coefs1, knots=knots, degree=degree, bsmooth=bsmooth, pcure=pcure1)
    h[,1] <- hsurvmspline(wane_period[1], alpha1, coefs=coefs1, knots=knots, degree=degree, bsmooth=bsmooth, pcure=pcure1)

    if (nints>1){ ## hazard interpolation needed
      ## weight for treated group in waned hazard.  w=1 is no waning, 0 is completely waned
      w <- (wane_period[2] - times) / (wane_period[2] - wane_period[1])
      for (i in 2:nints){
        logh1 <- hsurvmspline(times[i], alpha=alpha1, coefs=coefs1, pcure=pcure1,
                              knots=knots, degree=degree, bsmooth=bsmooth, log=TRUE)
        logh0 <- hsurvmspline(times[i], alpha=alpha0, coefs=coefs0, pcure=pcure0,
                              knots=knots, degree=degree, bsmooth=bsmooth, log=TRUE)
        h[,i] <- exp(w[i]*logh1 + (1 - w[i])*logh0)
        ## approximate cumulative hazard between time 0 and end of this interval
        H[,i] <- H[,i-1] + dt*h[,i-1]
      }
    }

    ## interpolate between the final cut point and x, assuming constant hazard in this interval
    mid <- (x > wane_period[1]) & (x < wane_period[2])
    w <- (wane_period[2] - times[ints[mid]]) / (wane_period[2] - wane_period[1])

    logh1 <- hsurvmspline(times[ints[mid]], alpha=alpha1[mid], coefs=coefs1[mid,,drop=FALSE], pcure=pcure1[mid],
                          knots=knots, degree=degree, bsmooth=bsmooth, log=TRUE)
    logh0 <- hsurvmspline(times[ints[mid]], alpha=alpha0[mid], coefs=coefs0[mid,,drop=FALSE], pcure=pcure0[mid],
                          knots=knots, degree=degree, bsmooth=bsmooth, log=TRUE)
    h_x <- exp(w*logh1 + (1-w)*logh0)

    hsurvmspline(0.1, alpha=alpha1[mid], coefs=coefs1[mid,,drop=FALSE], pcure=pcure1[mid], 
                 knots=knots, degree=degree, bsmooth=bsmooth, log=TRUE)
    dh_mid <- (x[mid] - times[ints[mid]]) * h_x
    ret[mid] <- H[cbind(which(mid), ints[mid])] + dh_mid

    ## exact cumulative hazard is known for intervals after the waning period
    late <- x >= wane_period[2]
    dh_late  <- Hsurvmspline(x[late], alpha=alpha0[late], coefs=coefs0[late,,drop=FALSE], pcure=pcure0[late], 
                             knots=knots, degree=degree, bsmooth=bsmooth) -
      Hsurvmspline(wane_period[2], alpha=alpha0[late], coefs=coefs0[late,,drop=FALSE], pcure=pcure0[late], 
                   knots=knots, degree=degree, bsmooth=bsmooth)

    ret[late] <- H[cbind(which(late), ints[late])] + dh_late
  }
  ret[offsetH>0] <- ret[offsetH>0] + offsetH[offsetH>0]

  ret
}

##' @rdname Survmspline_wane
##' @export
psurvmspline_wane <- function(q, alpha1, alpha0, coefs1, coefs0,
                              knots, degree=3, bsmooth=TRUE,
                              wane_period, wane_nt=10, lower.tail=TRUE,
                              pcure1=0, pcure0=0, offsetH=0, backhaz=NULL, log.p=FALSE){
  H <- Hsurvmspline_wane(x=q, alpha1=alpha1, alpha0=alpha0,
                         coefs1=coefs1, coefs0=coefs0,
                         knots=knots, degree=degree, bsmooth=bsmooth,
                         wane_period=wane_period, wane_nt=wane_nt,
                         pcure1=pcure1, pcure0=pcure0, offsetH=offsetH, backhaz=backhaz)
  ret <- if (lower.tail) 1 - exp(-H) else exp(-H)
  if (log.p) log(ret) else ret
}


##' @rdname Survmspline_wane
##' @export
hsurvmspline_wane <- function(x, alpha1, alpha0, coefs1, coefs0, knots, degree=3, bsmooth=TRUE,
                              wane_period, wane_nt=10, pcure1=0, pcure0=0, offseth=0, backhaz=NULL, log=FALSE){
  if (!is.null(backhaz)) offseth <- backhaz$hazard[findInterval(x, backhaz$time)]
  if (is.null(pcure1)) pcure1 <- 0
  if (is.null(pcure0)) pcure0 <- 0
  if (length(alpha0) != length(alpha1)) stop("lengths of `alpha0` and `alpha1` should be the same")
  att <- NULL
  d <- survmspline_dist_setup(q=x, alpha=alpha1, coefs=coefs1, knots=knots, pcure=pcure1, offseth=offseth)
  for (i in seq_along(d)) assign(names(d)[i], d[[i]])
  x <- q
  alpha1 <- d$alpha; coefs1 <- d$coefs; pcure1 <- d$pcure
  d0 <- survmspline_dist_setup(q=x, alpha=alpha0, coefs=coefs0, knots=knots, pcure=pcure0, offseth=offseth)
  alpha0 <- d0$alpha; coefs0 <- d0$coefs; pcure0 <- d0$pcure

  if (wane_period[1] >= wane_period[2]) stop("wane_period[2] should be greater than wane_period[1]")
  if (is.null(wane_nt)) wane_nt <- 10
  dt <- (wane_period[2] - wane_period[1]) / wane_nt

  ret <- numeric(length(x))
  early <- x <= wane_period[1]
  late <- x >= wane_period[2]
  ret[early] <- hsurvmspline(x[early], alpha=alpha1[early], coefs=coefs1[early,,drop=FALSE], pcure=pcure1[early], offseth=offseth[early],
                             knots=knots, degree=degree, bsmooth=bsmooth)
  ret[late] <- hsurvmspline(x[late],  alpha=alpha0[late], coefs=coefs0[late,,drop=FALSE], pcure=pcure0[late], offseth=offseth[late],
                            knots=knots, degree=degree, bsmooth=bsmooth)
  mid <- x > wane_period[1] & x < wane_period[2]
  if (any(mid)) {
    times <- c(0, seq(wane_period[1], wane_period[2], by=dt))
    ints <- findInterval(x[mid], times)
    w <- (wane_period[2] - times[ints]) / (wane_period[2] - wane_period[1])

    logh1 <- hsurvmspline(times[ints], alpha=alpha1[mid], coefs=coefs1[mid,,drop=FALSE], pcure=pcure1[mid], offseth=offseth[mid],
                 knots=knots, degree=degree, bsmooth=bsmooth, log=TRUE)
    logh0 <- hsurvmspline(times[ints], alpha=alpha0[mid], coefs=coefs0[mid,,drop=FALSE], pcure=pcure0[mid], offseth=offseth[mid],
                 knots=knots, degree=degree, bsmooth=bsmooth, log=TRUE)
    ret[mid] <- exp(w*logh1 + (1-w)*logh0)

  }

  ret[offseth>0] <- ret[offseth>0] + offseth[offseth>0]
  if (log) ret <- log(ret)

  attributes(ret) <- att
  ret
}

##' @rdname Survmspline_wane
##' @export
dsurvmspline_wane <- function(x, alpha1, alpha0, coefs1, coefs0, knots, degree=3, bsmooth=TRUE,
                              wane_period, wane_nt=10, pcure1=0, pcure0=0, offseth=0, offsetH=0, backhaz=NULL, log=FALSE){
  haz <- hsurvmspline_wane(x=x, alpha1=alpha1, alpha0=alpha0, coefs1=coefs1, coefs0=coefs0,
                           knots=knots, degree=degree, bsmooth=bsmooth,
                           wane_period=wane_period, wane_nt=wane_nt,
                           pcure1=pcure1, pcure0=pcure0, offseth=offseth,
                           backhaz=backhaz)
  surv <- psurvmspline_wane(q=x, alpha1=alpha1, alpha0=alpha0, coefs1=coefs1, coefs0=coefs0,
                            knots=knots, degree=degree, bsmooth=bsmooth,
                            wane_period=wane_period, wane_nt=wane_nt,
                            pcure1=pcure1, pcure0=pcure0, offsetH=offsetH,
                            backhaz=backhaz, lower.tail=FALSE)
  if (log) log(haz) + log(surv) else haz * surv
}


##' @rdname Survmspline_wane
##' @export
qsurvmspline_wane <- function(p, alpha1, alpha0, coefs1, coefs0, knots, degree=3, bsmooth=TRUE, lower.tail=TRUE, log.p=FALSE,
                              pcure1=0, pcure0=0, offsetH=0, backhaz=NULL, wane_period, wane_nt=10){
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (is.null(pcure1)) pcure1 <- 0
  if (is.null(pcure0)) pcure0 <- 0
  qgeneric(psurvmspline_wane, p=p, matargs=c("coefs1","coefs0"),
           scalarargs=c("knots","degree","bsmooth","backhaz","wane_period","wane_nt"),
           alpha1=alpha1, alpha0=alpha0, coefs1=coefs1, coefs0=coefs0,
           knots=knots, degree=degree, bsmooth=bsmooth,
           pcure1=pcure1, pcure0=pcure0,
           offsetH=offsetH, backhaz=backhaz, wane_period=wane_period, wane_nt=wane_nt)
}

##' @rdname Survmspline_wane
##' @export
rsurvmspline_wane <- function(n, alpha1, alpha0, coefs1, coefs0, knots, degree=3, bsmooth=TRUE,
                              pcure1=0, pcure0=0, offsetH=0, backhaz=NULL, wane_period, wane_nt=10){
  if (length(n) > 1) n <- length(n)
  qsurvmspline_wane(p=runif(n), alpha1=alpha1, alpha0=alpha0, coefs1=coefs1, coefs0=coefs0,
                    knots=knots, degree=degree, bsmooth=bsmooth,
                    pcure1=pcure1, pcure0=pcure0, offsetH=offsetH, backhaz=backhaz, wane_period=wane_period, wane_nt=wane_nt)
}

##' @rdname Survmspline_wane
##' @export
rmst_survmspline_wane = function(t, alpha1, alpha0, coefs1, coefs0,
                                 knots, degree=3, pcure1=0, pcure0=0,
                                 offsetH = 0, backhaz=NULL, bsmooth=TRUE,
                                 wane_period, wane_nt=10, disc_rate=0){
  if (is.null(pcure1)) pcure1 <- 0
  if (is.null(pcure0)) pcure0 <- 0
  rmst_generic(psurvmspline_wane, t, start=0,
               matargs = c("coefs1","coefs0"),
               unvectorised_args = c("knots","degree","backhaz","bsmooth","wane_period","wane_nt"),
               alpha1=alpha1, alpha0=alpha0, coefs1=coefs1, coefs0=coefs0, knots=knots, degree=degree,
               wane_period=wane_period, wane_nt=wane_nt, disc_rate=disc_rate,
               pcure1=pcure1, pcure0=pcure0,
               offsetH=offsetH, backhaz=backhaz, bsmooth=bsmooth)
}

##' @rdname Survmspline_wane
##' @export
mean_survmspline_wane = function(alpha1, alpha0, coefs1, coefs0, knots, degree=3, pcure1=0, pcure0=0,
                                 backhaz=NULL, bsmooth=TRUE,
                                 wane_period, wane_nt=10, disc_rate=0){
  if (is.null(pcure1)) pcure1 <- 0
  if (is.null(pcure0)) pcure0 <- 0
  nt <- if (is.matrix(coefs1)) nrow(coefs1) else 1
  rmst_generic(psurvmspline_wane, rep(Inf,nt), start=0,
               matargs = c("coefs1","coefs0"),
               unvectorised_args = c("knots","degree","backhaz","bsmooth","wane_period","wane_nt"),
               alpha1=alpha1, alpha0=alpha0, coefs1=coefs1, coefs0=coefs0, knots=knots, degree=degree,
               wane_period=wane_period, wane_nt=wane_nt, disc_rate=disc_rate,
               pcure1=pcure1, pcure0=pcure0, backhaz=backhaz, bsmooth=bsmooth)
}
