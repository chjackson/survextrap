##' Estimates of survival from a survextrap model
##'
##' @inheritParams print.survextrap
##' @inheritParams mean_survextrap
##'
##' @param times Vector of times at which to compute the estimates.
##'
##' @param tmax Maximum time at which to compute the estimates.  If
##' \code{times} is supplied, then this is ignored.  If \code{times} is not
##' supplied, then \code{times} is set to a set of 100 equally spaced
##' time points from 0 to \code{tmax}.  If both \code{tmax} and \code{times}
##' are not supplied, then \code{tmax} is set to the maximum follow up time
##' in the data.
##'
##' @param sample If \code{TRUE} then the MCMC samples are returned instead
##' of being summarised as a median and 95% credible intervals.
##'
##' @export
surv_survextrap <- function(x, newdata=NULL, times=NULL, tmax=NULL, niter=NULL, sample=FALSE) {
    if (is.null(newdata)) newdata <- x$mf_baseline
    pars <- get_pars(x, newdata=newdata, niter=niter)
    if (is.null(times)) times <- default_plottimes(x, tmax)
    nt <- length(times)
    niter <- attr(pars, "niter")
    nvals <- ncol(pars$linpreds)
    if (x$modelid==1){
        # times is [nt], alpha is [niter,1], coefs is [niter,K], linpreds is [niter,nvals]
        times_mat <- array(rep(times, niter*nvals), dim=c(nt, niter, nvals))
        linpreds_mat <- pars$linpreds[rep(1:niter, each=nt),]
        coefs_mat <- pars$coefs[rep(rep(1:niter, each=nt),nvals),]
        surv_sam <- psurvmspline(times_mat, linpreds_mat, coefs_mat,
                                  x$basehaz$knots, x$basehaz$degree, lower.tail=FALSE)
        dimnames(surv_sam) <- list(time=1:nt, iteration=1:niter, value=1:nvals)
    }
    else {
        ## TODO can we delete this
        times_mat <- matrix(times, ncol=length(times), nrow=niter, byrow=TRUE)
        surv_sam <- pweibull(times_mat,
                              shape=pars$coefs[,1], scale=exp(pars$alpha),
                              lower.tail=FALSE)
    }
    if (x$cure)  {
        cure_prob_mat <- array(rep(rep(pars$cure_prob, each=nt), nvals),
                               dim = c(nt, niter, nvals))
        surv_sam <- cure_prob_mat + (1 - cure_prob_mat)*surv_sam
    }
    if (!sample){
        surv <- apply(surv_sam, c(1,3), function(y)quantile(y, probs=c(0.5, 0.025, 0.975)))
        surv <- as.data.frame(matrix(surv, nrow=nt*nvals, ncol=3, byrow=TRUE))
        names(surv) <- c("median","lower","upper")
        surv <- cbind(times = rep(times, nvals), surv)
        if (!is.null(newdata))
            surv <- cbind(newdata[rep(1:nvals,each=nt),,drop=FALSE], surv)
    } else surv <- surv_sam
    attr(surv, "nvals") <- nvals
    surv
}

##' Estimates of hazard from a survextrap model
##'
##' @inheritParams print.survextrap
##' @inheritParams surv_survextrap
##'
##' @export
hazard_survextrap <- function(x, newdata=NULL, times=NULL, tmax=NULL, niter=NULL, sample=FALSE){
    if (is.null(newdata)) newdata <- x$mf_baseline
    pars <- get_pars(x, newdata=newdata, niter=niter)
    if (is.null(times)) times <- default_plottimes(x, tmax)
    nt <- length(times)
    niter <- attr(pars, "niter")
    nvals <- ncol(pars$linpreds)
    if (x$modelid==1){
        times_mat <- array(rep(times, niter*nvals), dim=c(nt, niter, nvals))
        linpreds_mat <- pars$linpreds[rep(1:niter, each=nt),]
        coefs_mat <- pars$coefs[rep(rep(1:niter, each=nt),nvals),]
        if (x$cure)
            logdens_sam <- dsurvmspline(times_mat, linpreds_mat, coefs_mat,
                                         x$basehaz$knots, x$basehaz$degree, log=TRUE)
        loghaz_sam <- hsurvmspline(times_mat, linpreds_mat, coefs_mat,
                                    x$basehaz$knots, x$basehaz$degree, log=TRUE)
        dimnames(loghaz_sam) <- list(time=1:nt, iteration=1:niter, value=1:nvals)
    }
    else  {
        ## TODO can we delete this
        times_mat <- matrix(times, ncol=length(times), nrow=niter, byrow=TRUE)
        logdens_sam <- stats::dweibull(times_mat,
                                 shape=pars$coefs[,1],
                                 scale=exp(pars$alpha),
                                 log=TRUE)
        loghaz_sam <- logdens_sam -
            stats::pweibull(times_mat, shape=pars$coefs[,1], scale=exp(pars$alpha),
                     lower.tail=FALSE, log.p=TRUE)
    }
    if (x$cure)  {
        cure_prob_mat <- array(rep(rep(pars$cure_prob, each=nt), nvals),
                               dim = c(nt, niter, nvals))
        loghaz_sam <- log(1 - cure_prob_mat) +  logdens_sam -
            log(surv_survextrap(x, times=times, sample=TRUE))
    }
    haz_sam <- exp(loghaz_sam)
    if (!sample){
        haz <- apply(haz_sam, c(1,3), function(y)quantile(y, probs=c(0.5, 0.025, 0.975)))
        haz <- as.data.frame(matrix(haz, nrow=nt*nvals, ncol=3, byrow=TRUE))
        names(haz) <- c("median","lower","upper")
        haz <- cbind(times = rep(times, nvals), haz)
        if (!is.null(newdata))
            haz <- cbind(newdata[rep(1:nvals,each=nt),,drop=FALSE], haz)
    } else haz <- haz_sam
    attr(haz, "nvals") <- nvals
    haz
}

##' Plot hazard curves from a survextrap model
##'
##' @inheritParams surv_survextrap
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
    haz <- hazard_survextrap(x, newdata=newdata, times=times, tmax=tmax, niter=niter)
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
##' @inheritParams surv_survextrap
##' @inheritParams plot_hazard
##'
##' @param km If \code{TRUE} then a Kaplan-Meier curve of the observed data is plotted,
##' using the results of \code{\link[survival:survfit]{survival::survfit()}} on the formula originally used
##' for the \code{survextrap} fit. 
##' By default, this is only done when there are no covariates or one factor covariate.
##'
##' @export
plot_survival <- function(x, newdata=NULL, times=NULL, tmax=NULL, km=NULL, niter=NULL,
                          ci=NULL, xlab="Time", ylab="Survival", line_size=1.5, ci_alpha=0.2){
    lower <- upper <- NULL
    surv <- surv_survextrap(x, newdata=newdata, times=times, tmax=tmax, niter=niter)
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
    if (is.null(tmax)) tmax <- max(c(x$eventtime, x$external$stop))
    times <- seq(tmin, tmax, by = (tmax - tmin) / nplot)
}

get_pars <- function(x, newdata=NULL, niter=NULL){
    stanmat <- as.matrix(x$stanfit)
    ## TODO error if there are no samples 
    if (is.null(niter)) niter <- nrow(stanmat)
    alpha    <- stanmat[1:niter, "alpha",  drop = FALSE]
    ms_coef_names <- sprintf("coefs[%s]",seq(x$nvars))
    coefs      <- stanmat[1:niter, ms_coef_names,  drop = FALSE]
    cure_prob <- if (x$cure) stanmat[1:niter, "cure_prob[1]"] else NULL
    linpreds <- get_linpreds(x=x, stanmat=stanmat, newdata=newdata, niter=niter)
    res <- nlist(alpha, linpreds, coefs, cure_prob)
    attr(res, "niter") <- niter
    res
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
