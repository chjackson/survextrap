# Based on the same functions in rstanarm but vastly simplified.
# Extensions to be added if there is a need

#' Prior distributions and options
#'
#' @name priors
#'
#' @description The functions described on this page are used to specify the prior distributions for the parameters in a survextrap model.
#'
#' @param location Prior location. For the normal distribution, this is the mean.
#' Defaults to 0
#'
#' @param scale Prior scale.  For the normal distribution, this is the standard deviation.
#' Defaults to 2.5.
#'
#' @param df Prior degrees of freedom (only for Student t distribution).
#'
#' @param shape1 First shape parameter (for Beta distribution, defaults to 1).
#'
#' @param shape2 Second shape parameter (for Beta distribution, defaults to 1).
#'
#' @param shape Shape parameter (for Gamma distribution, defaults to 2).
#'
#' @param rate Rate parameter (for Gamma distribution, defaults to 1).
#'
#' @seealso \code{\link{survextrap}}.
#'
#' @return A named list to be used internally by the \pkg{survextrap} model fitting functions.
#'
#'
NULL


#' @rdname priors
#' @export
p_normal <- function(location = 0, scale = 2.5) {
  validate_positive_parameter(scale)
  res <- nlist(dist = "normal", distid=1, location, scale, df=1)
  class(res) <- "prior"
  res
}

#' @rdname priors
#' @export
p_t <- function(location = 0, scale = 2.5, df = 1) {
  validate_positive_parameter(scale)
  validate_positive_parameter(df)
  res <- nlist(dist = "t", distid=2, location, scale, df)
  class(res) <- "prior"
  res
}

#' @rdname priors
#' @export
p_beta <- function(shape1 = 1, shape2 = 1){
    validate_positive_parameter(shape1)
    validate_positive_parameter(shape2)
    res <- nlist(dist = "beta", shape1, shape2)
    class(res) <- "prior"
    res
}

#' @rdname priors
#' @export
p_gamma <- function(shape = 2, rate = 1){
    validate_positive_parameter(shape)
    validate_positive_parameter(rate)
    res <- nlist(dist = "gamma", shape, rate)
    class(res) <- "prior"
    res
}


# internal ----------------------------------------------------------------

# Check for positive parameter (NULL is used in rstanarm)
#
# @param x The value to check.
# @return Either an error is thrown or \code{TRUE} is returned invisibly.
validate_positive_parameter <- function(x) {
  nm <- deparse(substitute(x))
  if (!is.null(x)) {
    if (!is.numeric(x))
      stop(nm, " should be NULL or numeric", call. = FALSE)
    if (any(x <= 0))
      stop(nm, " should be positive", call. = FALSE)
  }
  invisible(TRUE)
}



get_priors <- function(loghaz, loghr, smooth, cure, logor_cure, x, xcure, est_smooth,
                       nonprop, prior_sdnp, log_crude_event_rate){
  if (is.null(loghaz)) loghaz <- p_normal(log_crude_event_rate, 20)
  else validate_prior(loghaz)
  loghr <- get_prior_coveffs(loghr, x, "loghr") 
  smooth <- get_prior_smooth(smooth, est_smooth)
  validate_prior(cure)
  logor_cure <- get_prior_coveffs(logor_cure, xcure, "logor_cure")
  sdnp <- get_prior_sdnp(prior_sdnp, x, nonprop)
  nlist(loghaz, loghr, smooth, cure, logor_cure, sdnp)
}

validate_prior <- function(prior, priorname=NULL, element=NULL){
    if (is.null(priorname)) priorname <- deparse(substitute(prior))
    priorfullname <- sprintf("prior_%s", priorname)

    valid_priors <- list(loghaz = c("normal","t"),
                         loghr = c("normal","t"),
                         cure = "beta",
                         smooth = "gamma",
                         logor_cure = c("normal","t"),
                         sdnp = "gamma")

    element <- if (is.null(element)) "" else sprintf("[[%s]]",element)
    if (!inherits(prior, "prior")){
        stop(sprintf("`%s%s` should be a call to a prior constructor function such as %s(), see help(survextrap)",
                     priorfullname, element, valid_priors[[priorname]][1]))
    } else {
        if (!(prior$dist %in% valid_priors[[priorname]])){
            vp <- paste0("p_",valid_priors[[priorname]])
            vpstr <- if (length(vp)==1) vp else paste("one of", paste(vp,collapse=","))
            stop(sprintf("`%s%s` should be %s. p_%s was supplied", priorfullname, element, vpstr, prior$dist))
        }
    }
}

get_prior_smooth <- function(prior, est_smooth){
    validate_prior(prior, "smooth")
    if (!est_smooth) return(list(shape=aa(numeric()), rate=aa(numeric())))
    prior
}

## @param prior_list User-supplied prior specification
##
## Can either be like p_normal(0,1) or list(cov1=p_normal(0,2), cov2=p_normal(1, 3))
##
## @param x  List of information describing the linear model , as returned by make_x
get_prior_coveffs <- function(prior, x, modelname){
    priorname <- sprintf("prior_%s", deparse(substitute(prior)))
    if (x$ncovs==0)
        return(list(dist=aa(numeric()),distid=aa(numeric()),
                    location=aa(numeric()), scale=aa(numeric()),
                    df=aa(numeric())))
    else if (is.null(prior))
      prior <- p_normal(0, 2.5)
    prior_list <- validate_prior_bycov(prior, x, priorname)
    for (i in seq_len(x$ncovs)){
      validate_prior(prior_list[[i]], modelname, i)
    }
    dist <- sapply(prior_list, function(x)x$dist)
    distid <- sapply(prior_list, function(x)x$distid)
    location <- sapply(prior_list, function(x)x$location)
    scale <- sapply(prior_list, function(x)x$scale)
    df <- sapply(prior_list, function(x){if (is.null(x$df)) 1 else x$df})
    res <- data.frame(term=x$xnames, dist=dist, distid=distid,
                      location=location, scale=scale, df=df)
    rownames(res) <- NULL
    res
}

## Checks on priors specified as a list with one component per covariate
## If supplied as a prior for a single component, replicate that into a list 

validate_prior_bycov <- function(prior, x, priorname){
  if (inherits(prior, "prior")){
    prior_list <- rep(list(prior), x$ncovs)
    names(prior_list) <- x$xnames
  } else {
    if (!is.list(prior))
      stop(sprintf("%s should be a call (or list of calls) to a prior constructor function",
                   priorname))
    if (length(prior)!=x$ncovs){
      plural <- if (length(prior)>1) "s" else ""
      stop(sprintf("%s has %s component%s, but there are %s coefficients in the model",
                   priorname, length(prior), plural, x$ncovs))
    }
    if (!identical(sort(names(prior)), sort(x$xnames))){
      quoted_names <- sprintf("\"%s\"", x$xnames)
      stop(sprintf("names of %s do not match names of covariate coefficients: %s",
                   priorname, paste(quoted_names,collapse=",")))
    }
    prior_list <- prior
  }
  prior_list
}

get_prior_sdnp <- function(prior, x, nonprop){
  if (!nonprop) return(NULL)
  if (is.null(prior)) prior <- p_gamma(2, 1)
  prior_list <- validate_prior_bycov(prior, x, priorname="prior_sdnp")
  for (i in seq_len(x$ncovs)){
    validate_prior(prior_list[[i]], "sdnp", i)
  }
  shapes <- sapply(prior_list, function(x)x$shape)
  rates <- sapply(prior_list, function(x)x$rate)
  cbind(shapes, rates)
}
