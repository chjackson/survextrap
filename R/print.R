##' Print a fitted survextrap model
##'
##' @param x A fitted model object as returned by \code{\link{survextrap}}
##'
##' @param ... Other arguments (currently unused).
##'
##' @export
print.survextrap <- function(x, ...){
    cat("M-spline survival model\n")
    cat(sprintf("%s knots, degree %s, %s basis terms.\n", length(x$basehaz$knots), x$basehaz$degree, x$basehaz$nvars))
    cat(sprintf("Smoothness SD: %s\n", round(x$smooth_sd, 2)))
    print(x$stanfit)
}


##' Posterior summary statistics for parameters of survextrap models
##'
##' The intervals are 95\% credible intervals.
##'
##' Suggestions for what else to add to this are welcome.  Convergence diagnostics?
##'
##' @param object A fitted model object as returned by \code{\link{survextrap}}
##'
##' @param ... Other argument (currently unused).
##'
##' @export
summary.survextrap <- function(object, ...){
    i <- .value <- .variable <- variable <- term <- .lower <- .upper <- NULL
    sam <- get_draws(object)
    alpha <- tidybayes::gather_rvars(sam, alpha)
    coefs <- tidybayes::gather_rvars(sam, coefs[i])
    summ <- alpha %>% mutate(i=NA)
    if (object$ncovs>0){
        loghr <- tidybayes::gather_rvars(sam, loghr[i]) %>%
            mutate(term=object$xnames)
        hr <- loghr %>%
            mutate(.value = exp(.value),
                   .variable = "hr")
        summ <- summ %>%
            full_join(loghr, by=c(".variable",".value","i")) %>%
            full_join(hr, by=c(".variable",".value","i","term"))
    }
    summ <- summ %>%
        full_join(coefs, by=c(".variable",".value","i")) %>%
        mutate(sd = posterior::sd(.value)) %>%
        tidybayes::median_qi(.value) %>%
        rename(variable=.variable,
               median=.value) %>%
        select(variable, term, median,
               lower=.lower, upper=.upper, sd, i)
    if (object$modelid == "weibull"){
        ## TODO should change the parameter names to something more sensible if we
        ## keep the Weibull model in.  Shape, scale, consistently with dweibull
        summ <- summ %>%
            dplyr::filter(!(variable=="coefs" & i==2))
        if (object$ncovs==0) {
            summ <- summ %>% dplyr::select(-i)
            summ$variable[summ$variable=="loghr"] <- "logtaf"
            summ$variable[summ$variable=="hr"] <- "taf"
        }
    }
    summ
}

get_draws <- function(x){
    if (x$fit_method %in% c("mcmc","vb"))
        fit <- x$stanfit
    else if (x$fit_method=="opt")
        fit <- x$stanfit$theta_tilde
    else fit <- NULL
    posterior::as_draws_matrix(fit)
}

## TODO Other info rstanarm does in print or summary
## Model Info:
##  formula:         Surv(years, status) ~ 1
##  observations:    929
##  events:          407 (43.8%)
##  right censored:  522 (56.2%)
##  delayed entry:   no
##  algorithm:       sampling
##  sample:          4000 (posterior sample size)
##  priors:          see help('prior_summary')
