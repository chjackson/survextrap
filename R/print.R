## TODO indicate data used in print output, eg whether external

##' Print a fitted survextrap model
##'
##' @param x A fitted model object as returned by \code{\link{survextrap}}
##'
##' @param ... Other arguments (currently unused).
##'
##' This currently prints a statement of the fitted model, the prior distributions, and a table of summary statistics of the posterior distributions of the model parameters.  For descriptions of the parameters in this summary table, see \code{\link{summary.survextrap}}.
##'
##' @export
print.survextrap <- function(x, ...){
  cat("M-spline survival model\n")
  cat(sprintf("%s knots, degree %s, %s basis terms.\n", length(x$mspline$knots), x$mspline$degree, x$mspline$nvars))
  cat(sprintf("Smoothness SD: %s\n", if(x$est_smooth) "full Bayes" else round(x$smooth_sd, 2)))

  print_priors(x)
  cat("Posterior summary:\n")
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

##' Print the priors used in a fitted survextrap model
##'
##' @inheritParams print.survextrap
##'
##' @export
print_priors <- function(x){
    cat("Priors:\n")
    cat(sprintf("  Baseline log hazard scale: %s\n", form_prior_str(x$priors$loghaz)))
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


# List of parameters in survextrap models that summary statistics are calculated for
.parlist <- list(
  list(name="alpha", dimnames=NULL),
  list(name="coefs", dimnames="basis_num"),
  list(name="pcure", dimnames=NULL),
  list(name="loghr", dimnames="term"),
  list(name="hr", dimnames="term"),
  list(name="logor_cure", dimnames="term"),
  list(name="or_cure", dimnames="term"),
  list(name="smooth_sd", dimnames=NULL),
  list(name="nperr", dimnames=c("term","basis_num")),
  list(name="sd_np", dimnames="term")
)
# This is a data frame, whose second component ("dimnames") is in list column format
.parlist <- list2DF(list(name = sapply(.parlist, function(x)x[[1]]),
                         dimnames = sapply(.parlist, function(x)x[[2]]),
                         ndims = sapply(.parlist, function(x)length(x[[2]]))))

# Return a character vector with names of variables (including indices) in a Stan model
# These were stored in different places according to whether optimizing() or sampling() was used
get_parnames <- function(stanfit){
  if (inherits(stanfit,"stanfit"))
    names(stanfit)
  else names(stanfit$par)
}

# Get a regex pattern matching a Stan variable name of any dimension
# such as myname[1,3,2], myname[4,1], myname
# ndims is the number of dimensions
# The index into each dimension is captured to store separately
varname_pattern <- function(ndims){
  patt <- "^.+"
  if (ndims > 0)
    patt <- paste0(patt,
                   paste0("\\[",
                          paste(rep("(.+)", ndims), collapse=","), "\\]"))
  paste0(patt, "$")
}

# Build a data frame of parameter names and indices from a Stan model
# Indices with different meanings are stored in different variables
# Currently-used indices include
#  "basis_num" (for spline coefficient)
#  "term" (covariate described by a covariate effect)
parnames_to_df <- function(x){
  parname <- get_parnames(x$stanfit)
  pardf <- data.frame(parname)
  pardf$variable <- gsub("^(.+)\\[.+\\]$", "\\1", pardf$parname)
  udn <- unique(unlist(.parlist$dimnames))
  for (i in udn){
    pardf[[i]] <- NA
  }
  pardf <- pardf[pardf$variable %in% .parlist$name,]
  for (i in seq_len(nrow(.parlist))){
    thispar <- pardf$variable==.parlist$name[i]
    dn <- unlist(.parlist$dimnames[i])
    for (j in seq_along(dn)){
      index_name <- dn[j]
      pardf[[index_name]][thispar] <-
        as.numeric(
          gsub(varname_pattern(.parlist$ndims[i]), sprintf("\\%s",j),
               pardf$parname[thispar])
        )
    }
  }
  for (i in udn){
    if (all(is.na(pardf[[i]]))) pardf[[i]] <- NULL
  }
  has_xterm <- pardf$variable %in% c("loghr", "hr", "nperr", "sd_np")
  pardf$term[has_xterm] <- x$x$xnames[pardf$term[has_xterm]]
  pardf$term[pardf$variable %in% c("logor_cure","or_cure")] <- x$xcure$xnames
  pardf$basis_num[pardf$variable=="nperr"] <- pardf$basis_num[pardf$variable=="nperr"] + 1
  pardf$variable <- factor(pardf$variable, levels=.parlist$name)
  pardf <- pardf[order(pardf$variable),]
  pardf$variable <- as.character(pardf$variable)
  pardf
}

##' Posterior summary statistics for parameters of survextrap models
##'
##' By default, a median and 95% credible intervals are presented, alongside
##' the Rhat convergence diagnostic and the bulk effective sample size (see the \pkg{posterior} package).  For models
##' fitted by optimisation rather than MCMC, the posterior mode is also returned.
##'
##' Any other posterior summary can be computed if the appropriate function to compute it
##' is supplied here.
##'
##' @param object A fitted model object as returned by \code{\link{survextrap}}
##'
##' @param ... Other functions to compute a summary statistic from a vector
##' of posterior samples can be supplied here.   Many useful such functions are
##' provided with the \pkg{posterior} package. If the arguments are named,
##' then these are used to name the returned data frame.  See the examples below.
##'
##' @return A data frame (actually a tibble) of summary statistics for the model parameters is returned.
##'
##' The parameters, as indicated in the `variable` column, are:
##'
##' `alpha`: Average baseline log hazard.  If there are covariates, this describes the
##' average hazard with continuous
##' covariates set to zero, and factor covariates set
##' to their baseline levels.
##'
##' `coefs`: Coefficients of the M-spline basis terms.  If a non-proportional hazards model
##' was fitted, these are with covariates set to zero or baseline levels.
##'
##' `loghr`: Log hazard ratios for each covariate in the model. For cure
##' models, this refers to covariates on survival for uncured people.  For non-proportional hazards
##' models, these are the effects of covariates on the scale parameter `eta` and represent
##' "average" hazard ratios, from which departures are modelled (see the "methods"
##' vignette for a full description of this model).
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
##' `nperr`: Standardised departures from proportional hazards in the non-proportional hazards model, defined as \eqn{b^{(np)}_{ks} / \sigma^{(np)}_s} (see the "methods" vignette for definitions of these).
##'
##' `sd_np`: Smoothness standard deviations \eqn{\sigma^{(np)}_s} for the non-proportionality effects.
##'
##' @examples
##' mod <- survextrap(Surv(years, status) ~ rx, data=colons, fit_method="opt")
##' summary(mod)
##' summary(mod, mean)
##' summary(mod, mean, ess_tail=posterior::ess_tail)
##'
##' @export
summary.survextrap <- function(object, ...){
  summ <- parnames_to_df(object)
  if (object$fit_method=="opt")
    summ$mode <- object$stanfit$par[summ$parname]
  sam <- get_draws(object)
  sam <- sam[,summ$parname,drop=FALSE]

  mcmc_summ <- posterior::summarise_draws(sam,
                                          median,
                                          ~quantile(.x, probs=c(0.025, 0.975)),
                                          sd,
                                          rhat = posterior::rhat,
                                          ess_bulk = posterior::ess_bulk,
                                          ...)

  ## Remove strange classes ("pillar_num"  "pillar_vctr" "vctrs_vctr")
  for (i in which(sapply(mcmc_summ, is.numeric)))
    mcmc_summ[[i]] <- as.numeric(mcmc_summ[[i]])

  names(mcmc_summ)[names(mcmc_summ) == "variable"] <- "parname"
  names(mcmc_summ)[names(mcmc_summ) %in% c("2.5%","97.5%")] <- c("lower","upper")

  summ <- merge(summ, mcmc_summ, by="parname", sort=FALSE)
  summ$parname <- summ$index <- NULL

  tibble::as_tibble(summ)
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
    summ <- summary(object)
    summ <- summ[!(summ$variable %in% c("hr", "or", "coefs", "smooth_sd")), ]
    coefs <- summ$median
    term <- if ("term" %in% names(summ)) ifelse(is.na(summ$term), "", paste0("_", summ$term)) else ""
    names(coefs) <- paste0(summ$variable, term)
    coefs
}
