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
##' @param niter Number of MCMC iterations to use to compute credible
##' intervals.  Set to a low value to make this function quicker, at the cost of
##' some approximation error (which may not be important for plotting).
##'
##' @export
mean_survextrap <- function(x, newdata=NULL, niter=NULL){
    res <- rmst_survextrap(x, t=Inf, newdata=newdata, niter=niter)
    res$variable <- "mean"
    res$t <- NULL
    res
}

##' Restricted mean survival time
##'
##' @inheritParams mean_survextrap
##' @inheritParams surv_survextrap
##'
##' @param t Vector of times.  The restricted mean survival time up to each one of these times
##' will be computed.
##'
##' @export
rmst_survextrap <- function(x, t, newdata=NULL, niter=NULL){
    variable <- NULL
    if (is.null(newdata)) newdata <- x$mfbase
    pars <- get_pars(x, newdata=newdata, niter=niter)
    niter <- attr(pars, "niter")
    nt <- length(t)
    nvals <- ncol(pars$linpreds)
    res <- vector(nvals, mode="list")
    for (j in 1:nvals){
        resmat <- matrix(nrow=niter*nvals, ncol=nt)
        colnames(resmat) <- t
        for (i in 1:nt){
            resmat[,i] <- rmst_survmspline(t[i], pars$linpreds[,j], pars$coefs, x$basehaz$knots, x$basehaz$degree)
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

newdata_to_X <- function(newdata, x){
    form <- delete.response(terms(x$formula))
    if (is.null(newdata))
        X <- matrix(1, nrow=1, ncol=1)
    else
        X <- model.matrix(form, newdata, xlev=x$xlevs)
    X
}

get_linpreds <- function(x, stanmat, newdata=NULL, niter=NULL){
    if (x$ncovs > 0){
        X <- newdata_to_X(newdata, x) # nvals x npars
        X <- drop_intercept(X)
        X <- sweep(X, 2, x$xbar, FUN = "-")
    } else X <- NULL
    X <- cbind("(Intercept)"=1, X)
    if (is.null(niter)) niter <- nrow(stanmat)
    loghr_names <- grep("loghr", colnames(stanmat), value=TRUE)
    beta <- stanmat[1:niter, c("alpha", loghr_names)]
    beta %*% t(X) # niter x nvals
}

get_pars <- function(x, newdata=NULL, niter=NULL){
    stanmat <- get_draws(x)
    if (length(stanmat)==0) stop("Stan model does not contain samples")
    if (is.null(niter)) niter <- nrow(stanmat)
    if (is.null(newdata)) newdata <- x$mfbase
    alpha    <- stanmat[1:niter, "alpha",  drop = FALSE]
    gamma <- stanmat[1:niter, "gamma[1]", drop=FALSE]
    loghr_names <- sprintf("loghr[%s]",seq(x$ncovs))
    loghr <- if (x$ncovs>0) stanmat[1:niter, loghr_names, drop=FALSE] else NULL
    ms_coef_names <- sprintf("coefs[%s]",seq(x$nvars))
    coefs      <- stanmat[1:niter, ms_coef_names,  drop = FALSE]
    cure_prob <- if (x$cure) stanmat[1:niter, "cure_prob[1]"] else NULL
    linpreds <- get_linpreds(x=x, stanmat=stanmat, newdata=newdata, niter=niter)
    res <- nlist(alpha, gamma, loghr, linpreds, coefs, cure_prob)
    attr(res, "niter") <- niter
    res
}
