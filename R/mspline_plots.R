### Utilities for illustrating M-spline models


#' Get times and basis for an illustration of an M-spline with given knots.
#'
#' @param iknots  Internal knots.  If not supplied, \code{df} has to be specified, in which case
#' the default is \code{df - degree - 1} equally spaced knots between the boundary knots.
#'
#' @param bknots Boundary knots.   Defaults to `c(0,10)`.
#'
#' @param tmin Minimum plotting time.  Defaults to lower boundary knot.
#'
#' @param tmax Maximum plotting time.  Defaults to upper boundary knot.
#' 
#' @param degree Spline polynomial degree.
#'
#' @param df Desired number of basis terms, or "degrees of freedom" in the spline.
#' If \code{iknots} is not supplied, the number of internal knots is then chosen to satisfy this.
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
    if (is.null(iknots)) {
        nik <- df - degree  - 1
        iknots <- seq(bknots[1], bknots[2], length.out=nik+2)[-c(1,nik+2)]
    }
    validate_knots(iknots, name="iknots")
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
    } else if (length(coefs) != nvars) stop(sprintf("length of coefs is %s, should be %s", length(coefs), nvars))
    coefs <- coefs/sum(coefs)
    haz <- rowSums(basis * rep(coefs, each=nrow(basis))) * scale
    hazdf <- data.frame(haz=haz, time=time)
}


#' Plot a M-spline function, showing how it is built up from its basis
#'
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
#' @param plot If \code{TRUE} then a `ggplot2` plot object is returned.  If \code{FALSE} then the
#' underlying data are returned and no plot is done.
#'
#' @export
plot_mspline <- function(iknots=NULL, bknots=c(0,10), df=10, degree=3, coefs=NULL, scale=1, 
                         tmin=0, tmax=10, plot=TRUE){

    value <- term <- haz <- NULL
    s <- mspline_plotsetup(iknots=iknots, bknots=bknots, tmin=tmin, tmax=tmax, degree=degree, df=df)
    time <- s$time; basis <- s$basis

    ## TODO validate p length
    hazdf <- mspline_sum_basis(basis, coefs, scale, time)

    bdf <- as.data.frame(basis) %>%
        cbind(time=time) %>%
        pivot_longer(cols=all_of(1:ncol(basis)), names_to="term") %>%
        left_join(hazdf, by="time")

    if (plot)
        ggplot(bdf, aes(x=time, y=value, group=term)) +
            geom_line() +
            geom_line(data=hazdf, aes(x=time, y=haz), col="blue", inherit.aes = FALSE, lwd=1.5) +
            xlab("") +
            ylab("") +
            theme_minimal() +
            theme(panel.grid.minor = element_blank()) +
            scale_x_continuous(breaks=c(s$iknots, s$bknots)) +
            geom_vline(xintercept=c(s$iknots, s$bknots), col="blue", lwd=0.6, alpha=0.3)
    else list(basis=bdf, hazard=hazdf)
}


# Plot sample from prior distribution of hazards
#
# @inheritParams mspline_plotsetup
# @inheritParams plot_mspline_priorpred
#
#
mspline_priorpred_df <- function(iknots=NULL, bknots=c(0,10), df=10, degree=3,
                                 prior_mean, prior_sd=1, scale=1, scale_sd=1,
                                 tmin=0, tmax=10,
                                 nsim=10){
    ## TODO validate knots vs prior mean
    np <- length(prior_mean) + 1
    beta <- matrix(nrow=nsim, ncol=np)
    beta[,1] <- 0
    for (j in 2:np){
        beta[,j] <- rlogis(nsim, prior_mean[j-1], prior_sd)
    }
    coefs <- exp(beta) / rowSums(exp(beta))

    s <- mspline_plotsetup(iknots=iknots, bknots=bknots,
                           tmin=tmin, tmax=tmax, degree=degree, df=df)
    time <- s$time; basis <- s$basis; iknots <- s$iknots; bknots <- s$bknots
    hazlist <- vector(nsim, mode="list")
    scale <- exp(rnorm(nsim, log(scale), scale_sd))
    for (i in 1:nsim){
        hazlist[[i]] <- mspline_sum_basis(basis, coefs=coefs[i,], scale=scale[i], time=time)
        hazlist[[i]]$rep <- i
    }
    hazdf <- do.call("rbind", hazlist)
    hazdf$rep <- factor(hazdf$rep)
    attr(hazdf, "iknots") <- iknots
    attr(hazdf, "bknots") <- bknots
    hazdf
}

##' Generate a sample from the prior distribution of M-spline hazard curves implied by
##' a particular mean and variance for the baseline curve and scale parameter.
##'
##' See `vignette("methods")` for more information on this model. 
##'
##' @inheritParams mspline_plotsetup
##'
##' @param prior_mean Vector of means for logistic distributions.
##'
##' @param prior_sd Scale or SD parameter for logistic distributions (`scale` in `dllogis`)
##'
##' @param scale Prior mean on log scale for the log-normal prior for the scale parameter.
##'
##' @param scale_sd Prior SD on log scale for the log-normal prior for the scale parameter.
##'
##' @param nsim Number of simulations to draw
##'
##' @param plot If \code{TRUE} then a `ggplot2` plot object is returned.  If \code{FALSE} then the
##' underlying data are returned and no plot is done.
##'
##' @export
plot_mspline_priorpred <- function(iknots=NULL, bknots=c(0,10), df=10, degree=3,
                                   prior_mean, prior_sd=1, scale=1, scale_sd=1,
                                   tmin=0, tmax=10,
                                   nsim=10, plot=TRUE)
{
    haz <- NULL
    hazdf <- mspline_priorpred_df(iknots=iknots, bknots=bknots, df=df, degree=degree,
                                  prior_mean=prior_mean, prior_sd=prior_sd, scale=scale,
                                  scale_sd=scale_sd,
                                  tmin=tmin, tmax=tmax,
                                  nsim=nsim)
    iknots <- attr(hazdf, "iknots")
    bknots <- attr(hazdf, "bknots")
    if (plot)
    ggplot(hazdf, aes(x=time, y=haz, group=rep)) +
        geom_line(alpha=0.5) +
        scale_x_continuous(breaks=c(iknots, bknots)) +
        theme_minimal() +
        theme(panel.grid.minor = element_blank()) +
        geom_vline(xintercept = bknots, col="gray50") +
        xlab("Time") + ylab("Hazard")
    else
        hazdf
}

##' Estimate M-spline basis weights which give a constant function.
##'
##' This works by choosing the basis coefficients that minimise the
##' variance between log hazard values at different time points.
##' It is used in \code{\link{survextrap}} to choose the default prior mean
##' for the hazard function.
##' 
##' @param iknots Internal knots.
##'
##' @param bknots Boundary knots.
##'
##' @param times Times to use to construct the numerical calculation.
##' By default, this is 20 equally-spaced times between the boundary knots.
##'
##' @param degree Spline polynomial degree.
##'
##' @export
mspline_uniform_weights <- function(iknots, bknots, times=NULL, degree=3){
    if (is.null(times)) times <- seq(bknots[1], bknots[2], length.out=20)
    basis <- splines2::mSpline(times, knots = iknots, Boundary.knots = bknots,
                               degree = degree, intercept = TRUE)
    nvars <- ncol(basis)
    varloghaz <- function(logp){
        p <- exp(logp)
        haz <- rowSums(basis * rep(p, each=nrow(basis)))
        var(log(haz))
    }
    logp0 <- rep(0, nvars)
    opt <- optim(logp0, varloghaz, control=list(maxit=10000))
    res <- exp(opt$par)
    res / sum(res)
}

