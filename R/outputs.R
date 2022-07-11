##' Mean survival time
##'
##' @inheritParams print.survextrap
##'
##' @param newdata Data frame of covariate values to compute the output for.
##' If there are covariates in the model and this is not supplied,
##' the following default is used:
##'
##' (a) if the only covariate is one factor variable, then the output is computed
##' for each level of this factor.
##'
##' (b) if there are multiple covariates, or any numeric covariates, then the output
##' is computed at the mean of each numeric covariate in the original data, and at the
##' baseline level of each factor covariate.
##'
##' For cure models, then if there are covariates on the cure fraction and the uncured
##' survival, then the covariates on the cure fraction are used by default. 
##'
##' @param niter Number of MCMC iterations to use to compute credible
##' intervals.  Set to a low value to make this function quicker, at the cost of
##' some approximation error (which may not be important for plotting).
##'
##' @param ... Other options (currently unused).
##'
##' @export
mean.survextrap <- function(x, newdata=NULL, niter=NULL, ...){
    res <- rmst(x, t=Inf, newdata=newdata, niter=niter)
    res$variable <- "mean"
    res$t <- NULL
    res
}

default_newdata <- function(x, newdata){
    if (is.null(newdata)){
        if (x$ncurecovs > 0)
            newdata <- x$xcure$mfbase
        else if (x$ncovs > 0)
            newdata <- x$x$mfbase
    }
    newdata
}

##' Restricted mean survival time
##'
##' @inheritParams mean.survextrap
##'
##' @param t Vector of times.  The restricted mean survival time up to each one of these times
##' will be computed.
##'
##' @export
rmst <- function(x, t, newdata=NULL, niter=NULL){
    variable <- NULL
    newdata <- default_newdata(x, newdata)
    pars <- get_pars(x, newdata=newdata, niter=niter)
    niter <- attr(pars, "niter")
    nt <- length(t)
    nvals <- if (is.null(newdata)) 1 else nrow(newdata)
    res <- vector(nvals, mode="list")

    for (j in 1:nvals){
        resmat <- matrix(nrow=niter*nvals, ncol=nt)
        colnames(resmat) <- t
        for (i in 1:nt){
            if (x$modelid=="mspline")
                resmat[,i] <- rmst_survmspline(t[i], pars$alpha_user[,j], pars$coefs, x$basehaz$knots, x$basehaz$degree)
            else if (x$modelid=="weibull")
                resmat[,i] <- rmst_generic(pweibull, t[i], start=0, shape=pars$coefs[,1], scale=exp(pars$alpha_user[,j]))
            else stop("unknown modelid")
            if (x$cure){
                cureprob_mat <- matrix(rep(pars$cureprob_user, nt), nrow=niter*nvals, ncol=nt)
                tmat <- matrix(t, nrow=niter*nvals, ncol=nt, byrow=TRUE)
                resmat <- cureprob_mat*tmat + (1 - cureprob_mat)*resmat
            }
        }
        sample <- posterior::as_draws(resmat)
        res[[j]] <- summary(sample, median, ~quantile(.x, probs=c(0.025, 0.975))) %>%
            dplyr::rename(t = variable) %>%
            dplyr::mutate(t = as.numeric(t)) %>%
            dplyr::mutate(variable = "rmst") %>%
            dplyr::relocate(variable, .before=t)
    }
    res <- do.call("rbind", res)
    if (!is.null(newdata))
        res <- cbind(newdata[rep(1:nvals,each=nt),,drop=FALSE], res)
    attr(res, "sample") <- sample
    res
}

newdata_to_X <- function(newdata, x, formula=NULL, xlevs=NULL){
    if (is.null(xlevs)) xlevs <- x$x$xlevs
    if (is.null(formula)) formula <- x$formula
    form <- delete.response(terms(formula))
    if (is.null(newdata))
        X <- matrix(1, nrow=1, ncol=1)
    else
        X <- model.matrix(form, newdata, xlev=xlevs)
    X
}

get_alpha_bycovs <- function(x, stanmat, newdata=NULL, niter=NULL){
    if (is.null(newdata) && (x$ncovs > 0)) newdata <- x$x$mfbase
    if (NROW(newdata) > 0){
        X <- newdata_to_X(newdata, x, x$formula, x$x$xlevs) # nvals x npars
        X <- drop_intercept(X)
        X <- sweep(X, 2, x$x$xbar, FUN = "-")
    } else X <- NULL
    X <- cbind("(Intercept)"=1, X)
    if (is.null(niter)) niter <- nrow(stanmat)
    loghr_names <- grep("loghr", colnames(stanmat), value=TRUE)
    beta <- stanmat[1:niter, c("alpha", loghr_names)]
    beta %*% t(X) # niter x nvals
}

get_cureprob_bycovs <- function(x, stanmat, newdata=NULL, niter=NULL){
    if (is.null(newdata)) newdata <- x$xcure$mfbase
    pcure <- if (x$cure) stanmat[1:niter, "pcure[1]"] else NULL
    if (x$ncurecovs > 0){
        X <- newdata_to_X(newdata, x, x$cure_formula, x$xcure$xlevs) #
        X <- drop_intercept(X)
        X <- sweep(X, 2, x$xcure$xbar, FUN = "-")
        X <- cbind("(Intercept)"=1, X)
        if (is.null(niter)) niter <- nrow(stanmat)
        logor_names <- grep("logor_cure", colnames(stanmat), value=TRUE)
        logit_intercept <- qlogis(pcure)
        beta <- cbind(logit_intercept, stanmat[1:niter, c(logor_names)])
        pcure <- plogis(beta %*% t(X))
    }
    pcure
}

### TODO eta is inconsistently named in the code, log scale in stan, natural in vignette

get_pars <- function(x, newdata=NULL, niter=NULL){
    stanmat <- get_draws(x)
    if (length(stanmat)==0) stop("Stan model does not contain samples")
    if (is.null(niter)) niter <- nrow(stanmat)
    gamma <- stanmat[1:niter, "gamma[1]", drop=FALSE]
    loghr_names <- sprintf("loghr[%s]",seq(x$ncovs))
    loghr <- if (x$ncovs>0) stanmat[1:niter, loghr_names, drop=FALSE] else NULL
    ms_coef_names <- sprintf("coefs[%s]",seq(x$nvars))
    coefs      <- stanmat[1:niter, ms_coef_names,  drop = FALSE]

    ## cure probs for centered covariate values
    pcure <- if (x$cure) stanmat[1:niter, "pcure[1]"] else NULL
    ## cure probs for user-specified covariate values
    cureprob_user <- get_cureprob_bycovs(x=x, stanmat=stanmat, newdata=newdata, niter=niter)
    ## cure prob at covariate values of zero
    pcure0 <- get_cureprob_bycovs(x=x, stanmat=stanmat, newdata=x$xcure$mfzero, niter=niter)

    ## log intercept at centered covariate values (including centered factor contrasts)
    alpha    <- stanmat[1:niter, "alpha",  drop = FALSE]
    ## log intercept for user-specified covariate values
    alpha_user <- get_alpha_bycovs(x=x, stanmat=stanmat, newdata=newdata, niter=niter)
    ## log intercept at covariate values of zero
    alpha0 <- get_alpha_bycovs(x=x, stanmat=stanmat, newdata=x$x$mfzero, niter=niter)

    res <- nlist(alpha, alpha_user, alpha0, gamma, loghr, coefs, pcure, cureprob_user, pcure0)
    attr(res, "niter") <- niter
    res
}

##' Estimates of survival from a survextrap model
##'
##' @inheritParams print.survextrap
##' @inheritParams mean.survextrap
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
survival <- function(x, newdata=NULL, times=NULL, tmax=NULL, niter=NULL, sample=FALSE) {
    timedep_output(x, output="survival", newdata=newdata, times=times, tmax=tmax, niter=niter, sample=sample)
}

##' Estimates of hazard from a survextrap model
##'
##' @inheritParams print.survextrap
##' @inheritParams survival
##'
##' @export
hazard <- function(x, newdata=NULL, times=NULL, tmax=NULL, niter=NULL, sample=FALSE) {
    timedep_output(x, output="hazard", newdata=newdata, times=times, tmax=tmax, niter=niter, sample=sample)
}

timedep_output <- function(x, output="survival", newdata=NULL, times=NULL, tmax=NULL, niter=NULL, sample=FALSE){
    newdata <- default_newdata(x, newdata)
    pars <- get_pars(x, newdata=newdata, niter=niter)
    if (is.null(times)) times <- default_plottimes(x, tmax)
    nt <- length(times)
    niter <- attr(pars, "niter")
    nvals <- if (is.null(newdata)) 1 else nrow(newdata)
    times_arr <- array(rep(times, niter*nvals), dim=c(nt, niter, nvals))
    alpha_user_arr <- array(rep(pars$alpha_user, each=nt), dim=c(nt, niter, nvals))
    cureprob_arr <- if (x$cure) array(rep(pars$cureprob_user, each=nt), dim = c(nt, niter, nvals)) else NULL

    ## TODO cumulative hazard
    res_sam <- switch(output,
                      "survival" = survival_core(x, pars, times_arr, alpha_user_arr,
                                                 cureprob_arr, niter, nvals, nt),
                      "hazard" = hazard_core(x, pars, times_arr, alpha_user_arr,
                                             cureprob_arr, niter, nvals, nt)
                      )
    dimnames(res_sam) <- list(time=1:nt, iteration=1:niter, value=1:nvals)

    if (!sample){
        res <- apply(res_sam, c(1,3), function(y)quantile(y, probs=c(0.5, 0.025, 0.975)))
        res <- as.data.frame(matrix(res, nrow=nt*nvals, ncol=3, byrow=TRUE))
        names(res) <- c("median","lower","upper")
        res <- cbind(times = rep(times, nvals), res)
        if (!is.null(newdata))
            res <- cbind(newdata[rep(1:nvals,each=nt),,drop=FALSE], res)
    } else res <- res_sam
    attr(res, "nvals") <- nvals
    res
}

survival_core <- function(x, pars, times_arr, alpha_user_arr, cureprob_arr, niter, nvals, nt){
    if (x$modelid=="mspline"){
        coefs_mat <- pars$coefs[rep(rep(1:niter, each=nt),nvals),] # only for spline
        surv_sam <- psurvmspline(times_arr, alpha_user_arr, coefs_mat,
                                 x$basehaz$knots, x$basehaz$degree, lower.tail=FALSE)
    } else if (x$modelid=="weibull") {
        surv_sam <- pweibull(times_arr, shape=pars$coefs[,1],
                               scale=exp(alpha_user_arr), lower.tail=FALSE)
    }
    if (x$cure)  {
        surv_sam <- cureprob_arr + (1 - cureprob_arr)*surv_sam
    }
    surv_sam
}

hazard_core <- function(x, pars, times_arr, alpha_user_arr, cureprob_arr, niter, nvals, nt){
    if (x$modelid=="mspline"){
        coefs_arr <- pars$coefs[rep(rep(1:niter, each=nt),nvals),]
        if (x$cure)
            logdens_sam <- dsurvmspline(times_arr, alpha_user_arr, coefs_arr,
                                        x$basehaz$knots, x$basehaz$degree, log=TRUE)
        loghaz_sam <- hsurvmspline(times_arr, alpha_user_arr, coefs_arr,
                                   x$basehaz$knots, x$basehaz$degree, log=TRUE)
    } else if (x$modelid=="weibull") {
        logdens_sam <- stats::dweibull(times_arr, shape=pars$coefs[,1],
                                       scale=exp(alpha_user_arr), log=TRUE)
        loghaz_sam <- logdens_sam -
            stats::pweibull(times_arr, shape=pars$coefs[,1], scale=exp(pars$alpha),
                            lower.tail=FALSE, log.p=TRUE)
    }
    if (x$cure)  {
        loghaz_sam <- log(1 - cureprob_arr) +  logdens_sam -
            log(survival_core(x, pars, times_arr, alpha_user_arr, niter, nvals, nt))
    }
    haz_sam <- exp(loghaz_sam)
    haz_sam
}
