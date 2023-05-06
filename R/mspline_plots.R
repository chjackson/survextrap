##' Documentation for arguments used in illustration of M-splines
##'
##' @name mspline_args
##'
##' @param knots Vector of knot locations. If not supplied, \code{df} has to be specified, in which case
#' the default is \code{df - degree - 1} equally spaced knots between zero and the highest knot.
##' @param bknot Location of the highest knot
##'
##' @param degree Spline polynomial degree (defaults to 3)
##'
##' @param df Desired number of basis terms, or "degrees of freedom"
##'   in the spline.  If \code{knots} is not supplied, the number of
##'   knots is then chosen to satisfy this.
##'
##' @param bsmooth If \code{TRUE} then the function is constrained to
##'   also have zero derivative and second derivative at the boundary.
##'
NULL

mspline_even_knots <- function(knots=NULL, bknot, degree=3, df, bsmooth=TRUE){
  if (is.null(knots)) {
    nik <- if (bsmooth) df - degree  + 1 else df - degree  - 1
    knots <- seq(0, bknot, length.out=nik+2)[-1]
  }
  validate_knots(knots, name="knots")
  knots
}

#' Get basis for an illustration of an M-spline with given knots.
#'
#' @inheritParams mspline_args
#'
#' @param tmin Minimum plotting time.  Defaults to zero.
#'
#' @param tmax Maximum plotting time.  Defaults to the highest knot.
#'
mspline_plotsetup <- function(knots, bknot=10,
                              tmin=NULL, tmax=NULL,
                              degree=3, df=10, bsmooth=TRUE){
    if (is.null(tmin)) tmin <- 0
    if (is.null(tmax)) tmax <- bknot
    knots <- mspline_even_knots(knots, bknot, degree, df)
    time <- seq(tmin, tmax, length.out=1000)[-c(1,1000)]
    basis <- mspline_basis(times=time, knots=knots, degree=degree, bsmooth=bsmooth)
    basis
}

#' Plot a M-spline function, showing how it is built up from its basis
#'
#' @inheritParams mspline_plotdata
#'
#' @export
plot_mspline <- function(knots=NULL, bknot=10, df=10, degree=3, bsmooth=TRUE, coefs=NULL, scale=1,
                         tmin=0, tmax=10){
  haz <- term <- value <- NULL
    bdf <- mspline_plotdata(knots=knots, bknot=bknot, df=df, degree=degree, bsmooth=bsmooth, coefs=coefs,
                            scale=scale, tmin=tmin, tmax=tmax)
    ggplot(bdf, aes(x=time, y=value, group=term)) +
        geom_line() +
        geom_line(aes(x=time, y=haz), col="blue", inherit.aes = FALSE, lwd=1.5) +
        xlab("") +
        ylab("") +
        theme_minimal() +
        theme(panel.grid.minor = element_blank()) +
        scale_x_continuous(breaks=knots, limits=c(tmin, tmax)) +
        scale_y_continuous(limits=c(0, max(c(bdf$value, bdf$haz)))) +
        geom_vline(xintercept=knots, col="blue", lwd=0.6, alpha=0.3)
}

#' Data for plotting an M-spline function, showing how it is built up from its basis
#'
#' @inheritParams mspline_plotsetup
#'
#' @param coefs Coefficients of the spline basis terms.  These are normalised internally to sum to 1,
#' if they do not already sum to 1.
#'
#' @param scale Scale parameter. After computing the standard M-spline function as a weighted sum of the basis
#' terms, the function is multiplied by \code{scale}.   The log of the scale is the parameter called
#' \code{alpha} in the results of a `survextrap` model, the intercept of the linear model on the log hazard.
#'
#' @export
mspline_plotdata <- function(knots=NULL, bknot=10,
                             df=10, degree=3, bsmooth=TRUE, coefs=NULL, scale=1,
                             tmin=0, tmax=10){
    value <- term <- haz <- NULL
    basis <- mspline_plotsetup(knots=knots, bknot=bknot, tmin=tmin, tmax=tmax, degree=degree, df=df, bsmooth=bsmooth)
    time <- attr(basis, "times")
    hazdf <- data.frame(haz = mspline_sum_basis(basis, coefs, log(scale)),
                        time = time)
    bdf <- as.data.frame(basis)
    bdf$time <- time

    ## base R equivalent of pivot_longer(cols=all_of(1:ncol(basis)), names_to="term")
    bdf <- reshape(bdf, direction = "long",
                   varying = 1:ncol(basis),
                   v.names = "value",
                   timevar = "term")
    bdf <- bdf[order(bdf$id, bdf$term),]
    bdf$id <- NULL
    bdf <- merge(bdf, hazdf, by="time")

    attr(bdf,"knots") <- knots
    bdf
}


##' Generate and/or plot a sample from the prior distribution of M-spline hazard curves
##'
##' Generates and/or plots the hazard curves (as functions of time)
##' implied by a prior mean for the spline coefficients (a constant
##' hazard by default) and particular priors for the baseline log
##' hazard and smoothness standard deviation.
##'
##' @aliases mspline_priorpred plot_mspline_priorpred
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
##' @return \code{mspline_priorpred} returns a data frame of the
##'   samples, and \code{plot_mspline_priorpred} generates a plot.  No
##'   customisation options are provided for the plot function, which
##'   is just intended as a quick check.
##'
##' @name mspline_priorpred
##' @export
mspline_priorpred <- function(knots=NULL, df=10, degree=3, bsmooth=TRUE,
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

##' @rdname mspline_priorpred
##' @export
plot_mspline_priorpred <- function(knots=NULL, df=10, degree=3,
                                   coefs_mean = NULL,
                                   prior_hsd = p_gamma(2,1),
                                   prior_hscale = p_normal(0, 20),
                                   prior_loghr = NULL,
                                   X = NULL,
                                   prior_hrsd = p_gamma(2,1),
                                   tmin=0, tmax=NULL,
                                   nsim=10)
{
  haz <- NULL
  hazdf <- mspline_priorpred(knots=knots, df=df, degree=degree,
                             coefs_mean=coefs_mean,
                             prior_hsd=prior_hsd,
                             prior_hscale=prior_hscale,
                             prior_loghr=prior_loghr,
                             X = X,
                             prior_hrsd=prior_hrsd,
                             tmin=tmin, tmax=tmax, nsim=nsim)
  knots <- attr(hazdf, "knots")
  ggplot(hazdf, aes(x=time, y=haz, group=rep)) +
    geom_line(alpha=0.5) +
    scale_x_continuous(breaks=knots) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    geom_vline(xintercept = max(knots), col="gray50") +
    xlab("Time") + ylab("Hazard")
}
