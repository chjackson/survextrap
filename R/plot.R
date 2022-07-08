##' Plot hazard curves from a survextrap model
##'
##' @inheritParams survival
##'
##' @param ci If \code{TRUE} then credible intervals are drawn.  Defaults to
##' drawing the intervals if the plot shows the curve for only one covariate value.
##'
##' @param xlab X-axis label
##'
##' @param ylab Y-axis label
##'
##' @param line_size Passed to \code{\link{geom_line}}.
##'
##' @param ci_alpha Transparency for the credible interval ribbons
##'
##' @export
plot_hazard <- function(x, newdata=NULL, times=NULL, tmax=NULL, niter=NULL,
                        ci=NULL, xlab="Time", ylab="Hazard", line_size=1.5, ci_alpha=0.2){
    lower <- upper <- NULL # TODO do strings work
    haz <- hazard(x, newdata=newdata, times=times, tmax=tmax, niter=niter)
    knots <- x$basehaz$knots[x$basehaz$knots <= max(haz$times)]
    aes <- list(x="times", y="median")
    if (attr(haz, "nvals") > 1)
        aes <- c(aes, list(col = names(x$xlevs), group = names(x$xlevs)))
    geom_maps <- do.call("aes_string", aes)
    geom_ylab <- ggplot2::ylab(ylab)
    geom_xlab <- ggplot2::xlab(xlab)
    p <- ggplot(haz, mapping=geom_maps) +
        geom_ylab + geom_xlab +
        theme_minimal() +
        theme(panel.grid.minor = element_blank()) +
        geom_vline(xintercept=knots, col="blue", lwd=0.6, alpha=0.3) +
        geom_line(size=line_size)
    if (is.null(ci)) ci <- (attr(haz,"nvals")==1)
    if (ci)
        p <- p +
            geom_ribbon(aes(ymin=lower, ymax=upper), alpha=ci_alpha)
    p
}


##' Plot survival curves from a survextrap model
##'
##' @inheritParams survival
##' @inheritParams plot_hazard
##'
##' @param km If \code{TRUE} then a Kaplan-Meier curve of the observed data is plotted,
##' using the results of \code{\link[survival:survfit]{survival::survfit()}} on the formula originally used
##' for the \code{survextrap} fit.
##' By default, this is only done when there are no covariates or one factor covariate.
##'
##' The Kaplan-Meier estimates are returned in the \code{km} component of the fitted model object
##' returned by \code{\link{survextrap}}, if you want to construct your own plots like these.
##'
##'
##' @export
plot_survival <- function(x, newdata=NULL, times=NULL, tmax=NULL, km=NULL, niter=NULL,
                          ci=NULL, xlab="Time", ylab="Survival", line_size=1.5, ci_alpha=0.2){
    lower <- upper <- NULL
    surv <- survival(x, newdata=newdata, times=times, tmax=tmax, niter=niter)
    knots <- x$basehaz$knots[x$basehaz$knots <= max(surv$times)]
    aes <- list(x="times", y="median")
    if (attr(surv,"nvals") > 1)
        aes <- c(aes, list(col = names(x$xlevs), group = names(x$xlevs)))
    geom_maps <- do.call("aes_string", aes)
    geom_ylab <- ggplot2::ylab(ylab)
    geom_xlab <- ggplot2::xlab(xlab)
    g <- ggplot(surv, mapping=geom_maps) +
        ylim(0,1) +
        geom_ylab + geom_xlab +
        theme_minimal() +
        theme(panel.grid.minor = element_blank()) +
        geom_vline(xintercept=knots, col="blue", lwd=0.6, alpha=0.3) +
        geom_step(size=line_size)
    if (is.null(ci)) ci <- (attr(surv,"nvals")==1)
    if (ci)
        g <- g +
            geom_ribbon(aes(ymin=lower, ymax=upper), alpha=ci_alpha)

    if (is.null(km)) km <- one_factor_cov(x)
    if (km){
        aes <- list(x="time", y="surv")
        if (one_factor_cov(x)) {
            covname <- names(x$xlevs)
            aes <- c(aes, list(col=covname))
        }
        geom_maps <- do.call("aes_string", aes)
        g <- g + geom_step(data=x$km, mapping=geom_maps)
    }
    g
}

one_factor_cov <- function(x){
    (length(x$xinds$factor)<=1) && (length(x$xinds$numeric)==0)
}

default_plottimes <- function(x, tmax=NULL, nplot=100){
    tmin <- 0
    if (is.null(tmax)) tmax <- x$basehaz$bknots["upper"]
    times <- seq(tmin, tmax, by = (tmax - tmin) / nplot)
}

#' Plot method for survextrap model objects
#'
#' @param x  Fitted model object from \code{\link{survextrap}}.
#'
#' @param type `"survival"` for a plot of the survival function, `"hazard"` for the hazard function, against time.
#'
#' @param ... Additional arguments, passed on to \code{plot_hazard} and \code{plot_survival}.
#'
#' @import ggplot2
#'
#' @export
plot.survextrap <- function(x, type="hazsurv", ...){
    switch(type,
           "hazsurv" = plot_hazsurv(x, ...),
           "survival" = plot_survival(x, ...),
           "hazard" = plot_hazard(x, ...)
           )
}

plot_hazsurv <- function(x, ...){
    ps <- plot_survival(x, ...)
    ph <- plot_hazard(x, ...)
    gridExtra::grid.arrange(ps, ph, nrow=1)
}
