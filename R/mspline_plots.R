### Utilities for illustrating M-spline models


#' Get times and basis for an illustration of an M-spline with given knots.
#'
#' @inheritParams mspline_args
#'
#' @param tmin Minimum plotting time.  Defaults to lower boundary knot.
#'
#' @param tmax Maximum plotting time.  Defaults to upper boundary knot.
#'
#' @param extrap_model If `"constant"` then the hazard after the final boundary
#'  knot equals the hazard at this knot.  Otherwise the basis functions are simply
#'  extended outside boundary knot.
#'
mspline_plotsetup <- function(iknots, bknots=c(0,10),
                              tmin=NULL, tmax=NULL, degree=3, df=10,
                              extrap_model = "constant"){
    validate_knots(bknots, name="bknots")
    bknots <- sort(bknots)
    if (is.null(tmin)) tmin <- bknots[1]
    if (is.null(tmax)) tmax <- bknots[2]
    iknots <- mspline_default_iknots(iknots, bknots, degree, df)
    time <- seq(tmin, tmax, length.out=1000)[-c(1,length(time))]
    timeb <- time
    if (extrap_model=="constant"){
        timeb[timeb>bknots[2]] <- bknots[2]
        timeb[timeb<bknots[1]] <- bknots[1]
    }
    ## could also use basis_matrix for this, but allow alternative extrapolation method
    basis <- splines2::mSpline(timeb, knots = iknots, Boundary.knots = bknots,
                               degree = degree, intercept = TRUE)
    nlist(time, basis, iknots, bknots)
}

## Form the overall hazard function from summing basis terms
## Return as a tidy dataframe with time attached.
## ?? Do we need a version of this that returns a vector?
mspline_sum_basis <- function(basis, coefs=NULL, scale=1, time) {
    nvars <- ncol(basis)
    if (is.null(coefs)){
        coefs <- rep(1, nvars)
    } else if (length(coefs) != nvars) stop(sprintf("length of `coefs` is %s, should be %s", length(coefs), nvars))
    coefs <- coefs/sum(coefs)
    haz <- rowSums(basis * rep(coefs, each=nrow(basis))) * scale
    hazdf <- data.frame(haz=haz, time=time)
}


#' Plot a M-spline function, showing how it is built up from its basis
#'
#' @inheritParams mspline_plotdata
#'
#' @export
plot_mspline <- function(iknots=NULL, bknots=c(0,10), df=10, degree=3, coefs=NULL, scale=1,
                         tmin=0, tmax=10){
  haz <- term <- value <- NULL
    bdf <- mspline_plotdata(iknots=iknots, bknots=bknots, df=df, degree=degree, coefs=coefs,
                            scale=scale, tmin=tmin, tmax=tmax)
    knots <- c(attr(bdf,"iknots"), attr(bdf,"bknots"))
    ggplot(bdf, aes(x=time, y=value, group=term)) +
        geom_line() +
        geom_line(aes(x=time, y=haz), col="blue", inherit.aes = FALSE, lwd=1.5) +
        xlab("") +
        ylab("") +
        theme_minimal() +
        theme(panel.grid.minor = element_blank()) +
        scale_x_continuous(breaks=knots, limits=c(tmin, tmax)) +
        scale_y_continuous(breaks=knots, limits=c(0, max(c(bdf$value, bdf$haz)))) +
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
mspline_plotdata <- function(iknots=NULL, bknots=c(0,10), df=10, degree=3, coefs=NULL, scale=1,
                             tmin=0, tmax=10){
    value <- term <- haz <- NULL
    s <- mspline_plotsetup(iknots=iknots, bknots=bknots, tmin=tmin, tmax=tmax, degree=degree, df=df)
    time <- s$time; basis <- s$basis
    hazdf <- mspline_sum_basis(basis, coefs, scale, time)
    bdf <- as.data.frame(basis)
    bdf$time <- time

    ## base R equivalent of pivot_longer(cols=all_of(1:ncol(basis)), names_to="term")
    bdf <- reshape(bdf, direction = "long",
                   varying = list(as.character(1:ncol(basis))),
                   v.names = "value",
                   timevar = "term")
    bdf <- bdf[order(bdf$id, bdf$term),]
    bdf$id <- NULL
    bdf <- merge(bdf, hazdf, by="time")

    attr(bdf,"iknots") <- s$iknots
    attr(bdf,"bknots") <- s$bknots
    bdf
}


##' Generate and/or plot a sample from the prior distribution of M-spline hazard curves
##'
##' Generates and/or plots the hazard curves (as functions of time) implied by
##' a prior mean for the spline coefficients (a constant hazard by default)
##' and particular priors for the baseline log hazard and smoothness variance
##'
##' @aliases mspline_priorpred plot_mspline_priorpred
##'
##' @inheritParams survextrap
##' @inheritParams mspline_plotsetup
##'
##' @param nsim Number of simulations to draw
##'
##' @param nonprop Is the model a non-proportional hazards model (`TRUE` or `FALSE`).
##'
##' @param x Only required for models with covariates.  A list with components `ncovs` (number of covariates)
##' and `xnames` (names of these covariates).
##'
##' @param X Either a list, vector or one-row matrix with one covariate value.
##' Factors must be supplied as multiple numeric (0/1) contrasts. (TODO doc on how to build
##' the names - vignette)
##'
##' @param X0 An alternative covariate value (optional).  The hazard ratio between the hazards at X
##' and X0 will also be returned.
##'
##' @return \code{mspline_priorpred} returns a data frame of the samples, and
##' \code{plot_mspline_priorpred} generates a plot.
##'
##' @name mspline_priorpred
##' @export
mspline_priorpred <- function(iknots=NULL, bknots=c(0,10), df=10, degree=3,
                              coefs_mean = NULL,
                              prior_smooth = p_gamma(2,1),
                              prior_loghaz = NULL,
                              prior_loghr = NULL,
                              x = NULL,
                              X = NULL,
                              X0 = NULL, 
                              nonprop = FALSE,
                              prior_sdnp = p_gamma(2,1),
                              tmin=0, tmax=10,
                              nsim=10){
  ## TODO clean up args.  get knots from df, tmin and tmax.  How is this done?
  ## it is default_iknots.  what does that depend on?  just bknots and degree
  ## returns iknots.
  ##   if (is.null(mspline$iknots)) mspline$iknots <- mspline_default_iknots(...)
  
  s <- mspline_plotsetup(iknots=iknots, bknots=bknots,
                         tmin=tmin, tmax=tmax, degree=degree, df=df)
  time <- s$time; basis <- s$basis; iknots <- s$iknots; bknots <- s$bknots
  sam <- prior_sample(mspline = list(iknots=iknots, bknots=bknots, degree=degree),
                      coefs_mean=coefs_mean, prior_smooth=prior_smooth,
                      prior_loghaz=prior_loghaz,
                      prior_sdnp=prior_sdnp,
                      nonprop=nonprop,
                      x = x, X = X, nsim=nsim)
  hazlist <- vector(nsim, mode="list")
  for (i in 1:nsim){
    hazlist[[i]] <- mspline_sum_basis(basis, coefs=sam$coefs[i,],
                                      scale=exp(sam$loghaz[i]), time=time)
    hazlist[[i]]$rep <- i
  }
  hazdf <- do.call("rbind", hazlist)
  hazdf$mean <- 1/hazdf$haz
  hazdf$rep <- factor(hazdf$rep)

  if (!is.null(X0)){
    ## comparator hazard using the same parameter values but different covariate value
    sam <- prior_sample(mspline = list(iknots=iknots, bknots=bknots, degree=degree),
                        coefs_mean=coefs_mean, prior_smooth=prior_smooth,
                        prior_loghaz=prior_loghaz,
                        prior_sdnp=prior_sdnp,
                        nonprop=nonprop,
                        x = x, X = X0, nsim=nsim)
    hazlist <- vector(nsim, mode="list")
    for (i in 1:nsim){
      hazlist[[i]] <- mspline_sum_basis(basis, coefs=sam$coefs[i,],
                                        scale=exp(sam$loghaz[i]), time=time)
      hazlist[[i]]$rep <- i
    }
    hazdf0 <- do.call("rbind", hazlist)
    hazdf$haz0 <- hazdf0$haz
    hazdf$hr <- hazdf$haz / hazdf$haz0
  }

  attr(hazdf, "iknots") <- iknots
  attr(hazdf, "bknots") <- bknots
  hazdf
}

##' @rdname mspline_priorpred
##' @export
plot_mspline_priorpred <- function(iknots=NULL, bknots=c(0,10), df=10, degree=3,
                                   coefs_mean = NULL,
                                   prior_smooth = p_gamma(2,1),
                                   prior_loghaz = p_normal(0, 20),
                                   prior_loghr = NULL,
                                   x = NULL,
                                   X = NULL,
                                   nonprop = FALSE,
                                   prior_sdnp = p_gamma(2,1),
                                   tmin=0, tmax=10,
                                   nsim=10)
{
  haz <- NULL
  hazdf <- mspline_priorpred(iknots=iknots, bknots=bknots, df=df, degree=degree,
                             coefs_mean=coefs_mean,
                             prior_smooth=prior_smooth,
                             prior_loghaz=prior_loghaz,
                             tmin=0, tmax=10, nsim=10)
  iknots <- attr(hazdf, "iknots")
  bknots <- attr(hazdf, "bknots")
  ggplot(hazdf, aes(x=time, y=haz, group=rep)) +
    geom_line(alpha=0.5) +
    scale_x_continuous(breaks=c(iknots, bknots)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank()) +
    geom_vline(xintercept = bknots, col="gray50") +
    xlab("Time") + ylab("Hazard")
}
