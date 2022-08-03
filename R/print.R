## TODO indicate data used in print output, eg whether external 

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
    cat(sprintf("Smoothness SD: %s\n", if(x$est_smooth) "full Bayes" else round(x$smooth_sd, 2)))

    print_priors(x)
    cat("Posterior summary:\n")
## TODO convergence diagnostics
    print(summary(x))
}


form_prior_str <- function(prior){
    pars_db <- list(normal = c("location","scale"),
                    t = c("location","scale","df"),
                    beta = c("shape1","shape2"),
                    gamma = c("shape","rate"))
    pars <- prior[pars_db[[prior$dist]]]
    parstr <- paste(paste(names(pars), unlist(pars), sep="="), collapse=",")
    sprintf("%s(%s)", prior$dist, parstr)
}

print_priors <- function(x){
    cat("Priors:\n")
    cat(sprintf("  Baseline log hazard: %s\n", form_prior_str(x$priors$loghaz)))
    if (x$ncovs > 0){
        cat("  Log hazard ratios:\n")
        for (i in 1:x$ncovs){
            prior <- as.list(x$priors$loghr[i,])
            cat(sprintf("  %s: %s\n", prior$term, form_prior_str(prior)))
        }
    }
    if (x$cure)
        cat(sprintf("  Baseline cure probability: %s\n", form_prior_str(x$priors$cure)))
    if (x$ncurecovs > 0){
        cat("  Log odds ratios for cure probability:\n")
        for (i in 1:x$ncurecovs){
            prior <- as.list(x$priors$logor_cure[i,])
            cat(sprintf("  %s: %s\n", prior$term, form_prior_str(prior)))
        }
    }
    if (x$est_smooth){
        cat(sprintf("  Smoothness SD: %s\n", form_prior_str(x$priors$smooth)))
    }
}

##' Posterior summary statistics for parameters of survextrap models
##'
##' The intervals are 95\% credible intervals.
##'
##' Suggestions for what else to add to this are welcome.  Convergence diagnostics?
##' Mode for optim method.  Mean?
##' At least document how to get other arbitrary posterior summaries using
##' the posterior package.  Could even allow as an argument to this.
##'
##' @param object A fitted model object as returned by \code{\link{survextrap}}
##'
##' @param ... Other argument (currently unused).
##'
##' Document meanings of parameters
##'
##' `alpha`: Average baseline log hazard.  If there are covariates, this describes the
##' average hazard with continuous
##' covariates set to zero, and factor covariates set
##' to their baseline levels.
##'
##' `coefs`: Coefficients of the M-spline basis terms.
##'
##' `loghr`: Log hazard ratios for each covariate in the model. For cure
##' models, this refers to covariates on survival for uncured people.
##'
##' `hr`: Hazard ratios (the exponentials of `loghr`).
##'
##' `pcure`: Probability of cure (for cure models only).   If there are covariates
##' on cure, this parameter describes the probability of cure with continuous
##' covariates set to zero, and factor covariates set
##' to their baseline levels.
##'
##' `logor_cure`: Log odds ratio of cure for each covariate on cure.
##'
##' `or_cure`: Odds ratios of cure (the exponentials of `logor_cure`).
##'
##' @export
summary.survextrap <- function(object, ...){
    i <- .value <- .variable <- variable <- term <- .lower <- .upper <- logor_cure <- NULL
    sam <- get_draws(object)
    alpha <- tidybayes::gather_rvars(sam, alpha)
    coefs <- tidybayes::gather_rvars(sam, coefs[i])
    summ <- alpha %>% mutate(i=NA)
    if (object$cure){
        pcure <- tidybayes::gather_rvars(sam, pcure)
        summ <- summ %>% full_join(pcure, by=c(".variable",".value"))
    }
    if (object$ncovs>0){
        loghr <- tidybayes::gather_rvars(sam, loghr[i]) %>%
            mutate(term=object$x$xnames)
        hr <- loghr %>%
            mutate(.value = exp(.value),
                   .variable = "hr")
        summ <- summ %>%
            full_join(loghr, by=c(".variable",".value","i")) %>%
            full_join(hr, by=c(".variable",".value","i","term"))
    }
    if (object$ncurecovs>0){
        logor <- tidybayes::gather_rvars(sam, logor_cure[i]) %>%
            mutate(term=object$xcure$xnames)
        or <- logor %>%
            mutate(.value = exp(.value),
                   .variable = "or_cure")
        summ <- summ %>%
            full_join(logor, by=c(".variable",".value","i")) %>%
            full_join(or, by=c(".variable",".value","i","term"))
    }
    if (object$est_smooth){
        smooth_sd <- tidybayes::gather_rvars(sam, smooth_sd)
        summ <- summ %>% full_join(smooth_sd, by=c(".variable",".value"))
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

#' @export
coef.survextrap <- function(object, ...){
    variable <- NULL
    summ <- summary(object) %>%
        filter(!variable %in% c("hr", "or", "coefs", "smooth_sd"))
    coefs <- summ$median
    term <- if ("term" %in% names(summ)) ifelse(is.na(summ$term), "", paste0("_", summ$term)) else ""
    names(coefs) <- paste0(summ$variable, term)
    coefs
}
