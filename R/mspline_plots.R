
#' Get basis for an illustration of an M-spline with given knots.
#'
#' @inheritParams mspline_init
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
    mspline <- mspline_init(knots=knots, bknot=bknot,
                            degree=degree, df=df, bsmooth=bsmooth)
    time <- seq(tmin, tmax, length.out=1000)[-c(1,1000)]
    basis <- mspline_basis(times=time, knots=mspline$knots,
                           degree=mspline$degree, bsmooth=mspline$bsmooth)
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

##' @rdname prior_sample_hazard
##' @export
plot_prior_hazard <- function(knots=NULL, df=10, degree=3,
                              coefs_mean = NULL,
                              prior_hsd = p_gamma(2,1),
                              prior_hscale = p_normal(0, 20),
                              prior_loghr = NULL,
                              newdata = NULL,
                              formula = NULL,
                              prior_hrsd = p_gamma(2,1),
                              tmin=0, tmax=NULL,
                              nsim=10)
{
  haz <- NULL
  hazdf <- prior_sample_hazard(knots=knots, df=df, degree=degree,
                               coefs_mean=coefs_mean,
                               prior_hsd=prior_hsd,
                               prior_hscale=prior_hscale,
                               prior_loghr=prior_loghr,
                               newdata=newdata, formula=formula,
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
