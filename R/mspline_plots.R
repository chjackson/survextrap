### Utilities for illustrating M-spline models


#' Get times and basis for an illustration of an M-spline with given knots.
#'
#' @param knots  Internal knots.  If not supplied, df has to be specified.  Then defaults
#' to equally spaced knots between the boundary knots.
#'
#' @param bknots Boundary knots.   Defaults to `c(0,10)`.
#'
#' @param tmin Minimum plotting time.  Defaults to lower boundary knot
#'
#' @param tmax Maximum plotting time.  Defaults to upper boundary knot
#' 
#' @param degree Spline polynomial degree.
#'
#' @param df df TODO 
#'
#' @param extrap_model If `"constant"` then the hazard after the final boundary
#'  knot equals the hazard at this knot.  Otherwise the basis functions are simply
#'  extended outside boundary knot.
#'
mspline_plotsetup <- function(knots, bknots=c(0,10),
                              tmin=NULL, tmax=NULL, degree=3, df=10,
                              extrap_model = "constant"){
    if (is.null(tmin)) tmin <- bknots[1]
    if (is.null(tmax)) tmax <- bknots[2]
    if (is.null(knots)) {
        nik <- df - degree  - 1
        knots <- seq(bknots[1], bknots[2], length.out=nik+2)[-c(1,nik+2)]
    }
    if (is.null(tmax)) tmax <- max(bknots)
    time <- seq(tmin, tmax, length.out=1000)[-c(1,length(time))]
    timeb <- time
    if (extrap_model=="constant"){
        timeb[timeb>bknots[2]] <- bknots[2]
        timeb[timeb<bknots[1]] <- bknots[1]
    }
    ## could also use basis_matrix for this, but allow alternative extrapolation method 
    basis <- splines2::mSpline(timeb, knots = knots, Boundary.knots = bknots,
                               degree = degree, intercept = TRUE)
    nlist(time, basis, knots, bknots)
}

## Form the overall hazard function from summing basis terms
## Return as a tidy dataframe with time attached.
## ?? Do we need a version of this that returns a vector?
mspline_sum_basis <- function(basis, p=NULL, scale=1, time) {
    nvars <- ncol(basis)
    if (is.null(p)){
        p <- rep(1, nvars)
    } else if (length(p) != nvars) stop(sprintf("length of p is %s, should be %s", length(p), nvars))
    p <- p/sum(p)
    haz <- rowSums(basis * rep(p, each=nrow(basis))) * scale
    hazdf <- data.frame(haz=haz, time=time)
}

# Plot sample from prior distribution of hazards
#
# @param knots Vector of knots
#
# @param prior_mean Vector of means for logistic distributions.
#
# @param prior_sd Scale or SD parameter for logistic distributions (`scale` in `dllogis`)
#
mspline_priorpred_df <- function(knots=NULL, bknots=c(0,10), df=10, degree=3,
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
    p <- exp(beta) / rowSums(exp(beta))

    s <- mspline_plotsetup(knots=knots, bknots=bknots,
                           tmin=tmin, tmax=tmax, degree=degree, df=df)
    time <- s$time; basis <- s$basis; knots <- s$knots; bknots <- s$bknots
    hazlist <- vector(nsim, mode="list")
    scale <- exp(rnorm(nsim, log(scale), scale_sd))
    for (i in 1:nsim){
        hazlist[[i]] <- mspline_sum_basis(basis, p=p[i,], scale=scale[i], time=time)
        hazlist[[i]]$rep <- i
    }
    hazdf <- do.call("rbind", hazlist)
    hazdf$rep <- factor(hazdf$rep)
    attr(hazdf, "knots") <- knots
    attr(hazdf, "bknots") <- bknots
    hazdf
}

plot_mspline_priorpred <- function(knots=NULL, bknots=c(0,10), df=10, degree=3,
                                   prior_mean, prior_sd=1, scale=1, scale_sd=1,
                                   tmin=0, tmax=10,
                                   nsim=10)
{
    haz <- NULL
    hazdf <- mspline_priorpred_df(knots=knots, bknots=bknots, df=df, degree=degree,
                                  prior_mean=prior_mean, prior_sd=prior_sd, scale=scale,
                                  scale_sd=scale_sd,
                                  tmin=tmin, tmax=tmax,
                                  nsim=nsim)
    knots <- attr(hazdf, "knots")
    bknots <- attr(hazdf, "bknots")
    ggplot(hazdf, aes(x=time, y=haz, group=rep)) +
        geom_line(alpha=0.5) +
        scale_x_continuous(breaks=c(knots, bknots)) +
        theme_minimal() +
        theme(panel.grid.minor = element_blank()) +
        geom_vline(xintercept = bknots, col="gray50") +
        xlab("Time") + ylab("Hazard")
}

#' Plot a M-spline function, showing how it is built up from its basis
#'
#'
#' @param knots Internal knots
#'
#' @param bknots Boundary knots 
#'
#' @param p Weights TODO rename to coefs
#'
#' @param scale Scale parameter 
#'
#' @param degree Degree of polynomials defining the spline
#'
#' @param df TODO 
#'
#' @param tmin Minimum time to plot
#'
#' @param tmax Maximum time to plot to
#'
#' @param plot Plot 
#'
#' @export
plot_mspline <- function(knots=NULL, bknots=c(0,10), df=10, degree=3, p=NULL, scale=1, 
                         tmin=0, tmax=10, plot=TRUE){

    value <- term <- haz <- NULL
    s <- mspline_plotsetup(knots=knots, bknots=bknots, tmin=tmin, tmax=tmax, degree=degree, df=df)
    time <- s$time; basis <- s$basis

    ## TODO validate p length
    hazdf <- mspline_sum_basis(basis, p, scale, time)

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
            scale_x_continuous(breaks=c(knots, bknots)) +
            theme_minimal() +
            theme(panel.grid.minor = element_blank())
    else list(basis=bdf, hazard=hazdf)
}

##' Estimate M-spline basis weights which give a constant hazard
##'
##' @param knots TODO 
##'
##' @param bknots TODO 
##'
##' @param times TODO 
##'
##' @param degree TODO 
##'
##' @export
mspline_uniform_weights <- function(knots, bknots, times=NULL, degree=3){
    if (is.null(times)) times <- seq(bknots[1], bknots[2], length.out=20)
    basis <- splines2::mSpline(times, knots = knots, Boundary.knots = bknots,
                               degree = degree, intercept = TRUE)
    nvars <- ncol(basis)

    ## Function to calculate hazard given coefs.
    ## How to build in symmetry
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

