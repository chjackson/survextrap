## TODO more modularised plot options/defaults

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
##' @param line_size Passed to \code{\link{geom_line}}
##'
##' @param ci_alpha Transparency for the credible interval ribbons
##'
##' @param show_knots Show the locations of the spline knots as vertical lines
##'
##' @export
plot_hazard <- function(x, newdata=NULL, times=NULL, tmax=NULL, niter=NULL,
                        newdata0=NULL, wane_period=NULL, wane_nt=10,                        
                        ci=NULL, xlab="Time", ylab="Hazard",
                        line_size=1.5, ci_alpha=0.2, show_knots=FALSE){
    lower <- upper <- NULL
    newdata <- default_newdata(x, newdata)
    haz <- hazard(x, newdata=newdata, times=times, tmax=tmax, niter=niter,
                  newdata0=newdata0, wane_period=wane_period, wane_nt=wane_nt)
    knots <- x$mspline$knots[x$mspline$knots <= max(haz$times)]
    aes_list <- list(x=sym("times"), y=sym("median"))
    if (attr(haz, "nvals") > 1)
        aes_list <- c(aes_list, list(col = sym(names(newdata)), group = sym(names(newdata))))
    geom_ylab <- ggplot2::ylab(ylab)
    geom_xlab <- ggplot2::xlab(xlab)
    p <- ggplot(haz, mapping=aes(!!!aes_list)) + # requires ggplot2 3.4.0
        geom_ylab + geom_xlab +
        theme_minimal() +
        theme(panel.grid.minor = element_blank()) +
        geom_line(linewidth=line_size)
    if (show_knots) 
      p <- p + geom_vline(xintercept=knots, col="blue", linewidth=0.4*line_size, alpha=0.3)
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
                          newdata0=NULL, wane_period=NULL, wane_nt=10,                        
                          ci=NULL, xlab="Time", ylab="Survival",
                          line_size=1.5, ci_alpha=0.2, show_knots=FALSE){
    lower <- upper <- NULL
    if (is.null(km)) km <- one_factor_cov(x)
    newdata <- default_newdata(x, newdata)
    surv <- survival(x, newdata=newdata, times=times, tmax=tmax, niter=niter,
                     newdata0=newdata0, wane_period=wane_period, wane_nt=wane_nt)
    knots <- x$mspline$knots[x$mspline$knots <= max(surv$times)]
    aes_list <- list(x=sym("times"), y=sym("median"))
    if (attr(surv,"nvals") > 1)
        aes_list <- c(aes_list, list(col = sym(names(newdata)), group = sym(names(newdata))))
    geom_ylab <- ggplot2::ylab(ylab)
    geom_xlab <- ggplot2::xlab(xlab)

    ## requires ggplot2 3.4.0 (Nov 2022). though see 
    ## https://stackoverflow.com/questions/59578369/what-does-cant-use-at-top-level-mean-and-how-to-resolve-it
    ## for workarounds if this is too bleeding edge

    g <- ggplot(surv, mapping=aes(!!!aes_list)) + 
        ylim(0,1) +
        geom_ylab + geom_xlab +
        theme_minimal() +
        theme(panel.grid.minor = element_blank()) +
        geom_step(linewidth=line_size)
    if (show_knots)
      p <- p + geom_vline(xintercept=knots, col="blue", linewidth=0.4*line_size, alpha=0.3)
    if (is.null(ci)) ci <- (attr(surv,"nvals")==1)
    if (ci)
        g <- g +
            geom_ribbon(aes(ymin=lower, ymax=upper), alpha=ci_alpha)

    if (km){
        aes_list <- list(x=sym("time"), y=sym("surv"))
        if (one_factor_cov(x) && !is.null(newdata)) {
            aes_list <- c(aes_list, list(col=sym(names(newdata))))
        }
        g <- g + geom_step(data=x$km, mapping=aes(!!!aes_list))
    }
    g
}

one_factor_cov <- function(x){
    (length(x$x$factors)<=1) &&
        (length(x$x$numerics)==0) &&
        (x$ncurecovs == 0)
    ## if there are cure covs, make people choose for themselves what curves to draw
}

default_plottimes <- function(x, tmax=NULL, nplot=100){
    tmin <- 0
    if (is.null(tmax)) tmax <- x$mspline$bknots["upper"]
    times <- seq(tmin, tmax, by = (tmax - tmin) / nplot)
}

#' Plot method for survextrap model objects
#'
#' @inheritParams survival
#'
#' @param type `"survival"` for a plot of the survival function, `"hazard"` for the hazard function, against time.
#'
#' @param ... Additional arguments, passed on to \code{plot_hazard} and \code{plot_survival}.
#'
#' @import ggplot2
#'
#' @export
plot.survextrap <- function(x, type="hazsurv", newdata=NULL, ...){
    switch(type,
           "hazsurv" = plot_hazsurv(x, newdata=newdata, ...),
           "survival" = plot_survival(x, newdata=newdata, ...),
           "hazard" = plot_hazard(x, newdata=newdata, ...)
           )
}

plot_hazsurv <- function(x, newdata=NULL, ...){
  ps <- plot_survival(x, newdata=newdata, ...) +
    theme(legend.position = "none")
  ph <- plot_hazard(x, newdata=newdata, ...)
  ## Get the plot grids the same size, even when the legend is removed from the left one
  ps <- ggplotGrob(ps) 
  ph <- ggplotGrob(ph)
  g <- cbind(ps, ph)
  grid::grid.newpage()
  grid::grid.draw(g)
}

##' Plot hazard ratio against time from a survextrap model
##'
##' Intended for use with non-proportional hazards models (\code{survextrap(...,nonprop=TRUE)}).
##'
##'
##' @inheritParams hazard_ratio
##' @inheritParams plot_hazard
##'
##' @export
plot_hazard_ratio <- function(x, newdata=NULL, times=NULL, tmax=NULL, niter=NULL,
                              ci=TRUE, xlab="Time", ylab="Hazard ratio",
                              line_size=1.5, ci_alpha=0.2){
    lower <- upper <- NULL
    hr <- hazard_ratio(x, newdata=newdata, times=times, tmax=tmax, niter=niter)
    knots <- x$mspline$knots[x$mspline$knots <= max(hr$times)]
    aes_list <- list(x=sym("times"), y=sym("median"))
    geom_ylab <- ggplot2::ylab(ylab)
    geom_xlab <- ggplot2::xlab(xlab)
    p <- ggplot(hr, mapping=aes(!!!aes_list)) +
        geom_ylab + geom_xlab +
        theme_minimal() +
        theme(panel.grid.minor = element_blank()) +
        geom_vline(xintercept=knots, col="blue", lwd=0.4*line_size, alpha=0.3) +
        geom_hline(yintercept=1, col="gray30") +
        geom_line(linewidth=line_size)
    if (ci)
        p <- p +
            geom_ribbon(aes(ymin=lower, ymax=upper), alpha=ci_alpha)
    p
}
