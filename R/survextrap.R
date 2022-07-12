## TODO
## use capital X for covariates not x
## Document / error that variables should be in data frame not in working env

#' survextrap
#'
#' @param formula  A survival formula in standard R formula syntax.  Covariates included
#' on the right hand side of the formula with be modelled with proportional hazards.
#'
#' @param data Data frame containing variables in `formula`.
#'
#' @param external External data as a data frame of aggregate survival counts with columns:
#'
#' `start`: Start time
#'
#' `start`: Follow-up time
#'
#' `n`: Number of people alive at `start`
#'
#' `r`: Number of those people who are still alive at `stop`
#'
#' It is intended to add facilities to the package to produce this format of data from
#' other common forms of external data, e.g. registry or population data, or elicited judgements
#' about conditional survival
#'
#' @param smooth_sd Smoothing parameter estimation.
#'
#' `"bayes"`: the smoothing parameter is estimated by full Bayes.
#'
#'  `"eb"`: empirical Bayes is used.
#'
#' Alternatively, if a number is supplied here, then the smoothing parameter is fixed to this number.
#'
#' @param coefs_mean Spline basis coefficients that define the prior mean for the hazard function. By
#' default, these are set to values that define a constant hazard function.
#'
#' @param cure If `TRUE`, a mixture cure model is used.
#'
#' @param cure_prior Beta shape parameters for the cure probability (vector of two).
#'
#' @param basehaz_ops A list of control parameters defining the spline model.
#'
#' `iknots`: Internal knots.  If this is not supplied, then the number of knots is taken from `df`,
#' and their location is taken from equally-spaced quantiles of the observed event times (concatenated
#' with the distinct follow-up times in the external data).
#'
#' `bknots`: Boundary knots.  If this is not supplied, the boundary knots are set to zero and the
#' maximum event time (including the follow-up times in the external data).
#'
#' `df`: Degrees of freedom, i.e. the number of parameters (or basis terms) that define the hazard
#' as a spline function of time.  Defaults to 10.  Note this does not necessarily overfit because
#' the function is smoothed through the prior.
#'
#' `degree`: Polynomial degree used for the basis function. The default is 3, giving a cubic.
#'
#' @param modelid \code{"mspline"} for the default M-spline model.
#'
#' Only current alternative is \code{"weibull"} for a Weibull accelerated failure time model.
#' Just included for the purpose of package testing. Not fully tested, and not recommended for use in
#' practice for survival extrapolation.
#'
#' @param fit_method Method from \pkg{rstan} used to fit the model.  Defaults to MCMC.
#'
#'  If \code{fit_method="mcmc"} then a sample from the posterior is drawn using Markov Chain Monte Carlo
##' sampling, via \pkg{rstan}'s \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}} function.
##' This is the default.  It is the most accurate but the slowest method.
##'
##'   If \code{fit_method="opt"}, then instead of an MCMC sample from the posterior,
##'   `disbayes` returns the posterior mode calculated using optimisation, via
##'   \pkg{rstan}'s \code{\link[rstan:stanmodel-method-optimizing]{rstan::optimizing()}} function.
##'   A sample from a normal approximation to the (real-line-transformed)
##'   posterior distribution is drawn in order to obtain credible intervals.
##'
##'   If \code{fit_method="vb"}, then variational Bayes methods are used, via \pkg{rstan}'s
##'   \code{\link[rstan:stanmodel-method-vb]{rstan::vb()}} function.  This is labelled as "experimental" by
##'   \pkg{rstan}.  It might give a better approximation to the posterior
##'   than \code{fit_method="opt"}, but has not been investigated much for `disbayes` models.
##'
##' @param loo Compute leave-one-out cross-validation statistics.  This is done by default for MCMC fits,
##' and not for optimisation or VB fits. Set to \code{FALSE} to save a bit of run time.
##' If these statistics are computed, then they are returned in the \code{loo} component of the
##' object returned by \code{survextrap}.
##'
#' @param ... Additional arguments to supply to control the Stan fit, passed to the appropriate
#' \pkg{rstan} function.
#'
#'
#' @export
survextrap <- function(formula,
                       data,
                       external=NULL,
                       smooth_sd = "bayes",
                       coefs_mean = NULL,
                       cure = FALSE,
                       cure_prior = c(1,1),
                       basehaz_ops = NULL,
                       modelid = "mspline",
                       fit_method = "mcmc",
                       loo = (fit_method=="mcmc"),
                       ...)
{

    Terms <- terms(formula) # no random effects for now
    mf <- stats::model.frame(Terms, data)

    t_beg <- make_t(mf, type = "beg") # entry time
    t_end <- make_t(mf, type = "end") # exit  time
    t_upp <- make_t(mf, type = "upp") # upper time for interval censoring [ not implemented yet ]
    status <- make_d(mf)
    t_event <- aa(t_end[status == 1]) # exact event time
    t_rcens <- aa(t_end[status == 0]) # right censoring time
    nevent <- sum(status == 1)
    nrcens <- sum(status == 0)
    ind_event <- which(status==1)
    ind_rcens <- which(status==0)

    x <- make_x(formula, mf)
    ncovs <- x$ncovs
    x_event <- x$x_centered[ind_event, , drop = FALSE]
    x_rcens <- x$x_centered[ind_rcens, , drop = FALSE]

    xcure <- make_xcure(cure,data=data)
    cure <- xcure$cure
    cure_formula <- xcure$cure_formula
    ncurecovs <- xcure$ncovs
    xcure_event <- xcure$x_centered[ind_event, , drop = FALSE]
    xcure_rcens <- xcure$x_centered[ind_rcens, , drop = FALSE]

    t_tmp <- sum(rowMeans(cbind(t_end, t_upp), na.rm = TRUE) - t_beg)
    d_tmp <- sum(!status == 0)
    log_crude_event_rate <- log(d_tmp / t_tmp)
    if (is.infinite(log_crude_event_rate))
        log_crude_event_rate <- 0 # avoids error when there are zero events

    external <- parse_external(external, formula, ncovs, xlevs=x$xlevs, xbar=x$xbar,
                               cure_formula, ncurecovs, xcure$xlevs, xcure$xbar)
    t_ext_stop <- aa(external$stop)
    t_ext_start <- aa(external$start)
    r_ext <- aa(external$r)
    n_ext <- aa(external$n)
    nextern <- external$nextern
    x_ext <- external$x_centered
    xcure_ext <- external$xcure_centered
    tmax <- max(c(t_end,t_upp,external$tmax), na.rm = TRUE)

    basehaz <- make_basehaz(basehaz_ops    = basehaz_ops,
                            times          = t_end,
                            times_ext      = unique(c(t_ext_start, t_ext_stop)),
                            status         = status,
                            tmin          = min(t_beg),
                            tmax          = tmax)

    basis_event  <- make_basis(t_event, basehaz)
    ibasis_event <- make_basis(t_event, basehaz, integrate = TRUE)
    ibasis_rcens <- make_basis(t_rcens, basehaz, integrate = TRUE)
    nvars <- basehaz$nvars

    if (is.null(coefs_mean)){
        coefs_mean <- mspline_uniform_weights(iknots = basehaz$iknots, bknots=basehaz$bknots)
    } else {
        ## TODO validate
    }
    beta_mean <- log(coefs_mean[-1] / coefs_mean[1])
    modelids <- c("mspline", "weibull")
    if (!(modelid %in% modelids)) stop("modelid should be mspline or weibull")
    modelid_num <- match(modelid, modelids)

    est_smooth <- (smooth_sd == "bayes")
    if (est_smooth) smooth_sd <- 1

    if (modelid=="weibull"){
        ## Weibull model.
        nvars <- 2 # Number of free parameters, consistent with M-spline.  M-spline model has an extra "1 minus sum of rest" parameter
        basis_event <- array(t_event, dim=c(length(t_event), 2)) # just first col of these used.
        ibasis_event <- array(t_event, dim=c(length(t_event), 2))
        ibasis_rcens <- array(t_rcens, dim=c(length(t_rcens), 2))
        beta_mean <- aa(1) # prior mean for log shape parameter. shape is coefs[1], log scale is eta[1]
        est_smooth <- FALSE
    }
    ibasis_ext_stop <- if (nextern>0) make_basis(t_ext_stop, basehaz, integrate = TRUE) else matrix(nrow=0, ncol=nvars)
    ibasis_ext_start <- if (nextern>0) make_basis(t_ext_start, basehaz, integrate = TRUE) else matrix(nrow=0, ncol=nvars)

    standata <- nlist(nevent, nrcens, nvars, nextern, ncovs,
                      log_crude_event_rate,
                      basis_event, ibasis_event, ibasis_rcens,
                      ibasis_ext_stop, ibasis_ext_start,
                      x_event, x_rcens,
                      r_ext, n_ext, x_ext,
                      ncurecovs, xcure_event, xcure_rcens, xcure_ext,
                      beta_mean,
                      est_smooth,
                      cure,
                      cure_shape = cure_prior, ## TODO standardise name, mean+ESS param
                      modelid = modelid_num
                      )
    pcure_init <- if (cure) 0.5 else numeric()
    staninit <- list(gamma = aa(0),
                     loghr = aa(rep(0, standata$ncovs)),
                     beta_err = aa(rep(0, standata$nvars)),
                     smooth_sd = aa(if(standata$est_smooth) smooth_sd else numeric()),
                     pcure = aa(pcure_init))
    if (identical(smooth_sd, "eb")){
        smooth_sd <- eb_smoothness(standata, staninit)
    }
    standata$smooth_sd_fixed <- if (est_smooth) aa(numeric()) else aa(smooth_sd)

    if (fit_method=="opt")
        fits <- rstan::optimizing(stanmodels$survextrap, data=standata, init=staninit,
                                  hessian=TRUE, draws=1000) # TODO defaults
    else if (fit_method=="mcmc")
        fits <- rstan::sampling(stanmodels$survextrap, data=standata,
                                pars = "beta", include=FALSE,
                                ...)
    else if (fit_method=="vb")
        fits <- rstan::vb(stanmodels$survextrap, data=standata,
                          pars = "beta", include=FALSE,
                          ...)
    else stop(sprintf("Unknown fit_method: %s",fit_method))

    km <- survminer::surv_summary(survfit(formula, data=data), data=data)

    misc_keep <- nlist(formula, stanfit=fits, fit_method, cure_formula)
    standata_keep <- standata[c("nvars","ncovs","log_crude_event_rate","ncurecovs")]
    model_keep <- nlist(cure, modelid)
    spline_keep <- nlist(basehaz)
    covinfo_names <- c("xnames","xlevs","xinds","xbar","mfbase","mfzero")
    x <- list(x = x[covinfo_names])
    xcure <- list(xcure = xcure[covinfo_names])
    prior_keep <- nlist(coefs_mean, smooth_sd)
    res <- c(misc_keep, standata_keep, model_keep, spline_keep, x, xcure, prior_keep, nlist(km))

    class(res) <- "survextrap"
    if (loo) res$loo <- loo_survextrap(res, standata)
    res
}


make_x <- function(formula, mf, xlevs=NULL){
    x <- model.matrix(formula, mf, xlev = xlevs)
    x <- drop_intercept(x)
    ncovs <- NCOL(x)
    xbar <- aa(colMeans(x))
    x_centered <- sweep(x, 2, xbar, FUN = "-")
    xnames <- colnames(x)

    has_response <- attr(terms(mf), "response") != 0
    if (has_response)
        mf <- mf[-attr(terms(mf), "response")]
    if (ncol(mf) > 0){
        xinds <- list(
            factor = which(sapply(mf, is.factor)),
            numeric = which(sapply(mf, is.numeric))
        )
        xlevs <- lapply(mf[xinds$factor], levels)
        if ((length(xinds$factor)==1) && (length(xinds$numeric)==0)){
            base_levs <- xlevs
        }
        else if (length(xinds$factor) > 0)
            base_levs <- lapply(xlevs, function(x)x[1])
        else base_levs <- NULL
        numeric_means <- lapply(mf[xinds$numeric], mean)
        numeric_zeros <- lapply(mf[xinds$numeric], function(x)0)
        mfbase <- as.data.frame(c(numeric_means, base_levs))
        mfzero <- as.data.frame(c(numeric_zeros, base_levs))
    } else xlevs <- xinds <- mfbase <- mfzero <- NULL

    nlist(x, x_centered, xbar, N = NROW(x),
          ncovs = NCOL(x), xnames,
          xlevs, xinds, mfbase, mfzero)
}


eb_smoothness <- function(standata, staninit){
    standata$est_smooth <- 1
    standata$smooth_sd_fixed  <- aa(numeric()) # dummy
    fits <- rstan::optimizing(stanmodels$survextrap, data=standata, init=staninit, hessian=FALSE, verbose=TRUE)
    if (fits$return_code==0){
        smooth_sd <- fits$par["smooth_sd[1]"]
    } else {
        warning("Empirical Bayes estimation of smoothness parameter failed, continuing with default")
        smooth_sd <- 1
    }
    smooth_sd
}

validate_formula <- function(formula, needs_response = TRUE) {

  if (!inherits(formula, "formula")) {
    stop2("'formula' must be a formula.")
  }

  if (needs_response) {
    len <- length(formula)
    if (len < 3) {
      stop2("'formula' must contain a response.")
    }
  }
  as.formula(formula)
}

# Check object is a Surv object with a valid type
#
# @param x A Surv object. That is, the LHS of a formula as evaluated in a
#   data frame environment.
# @param ok_types A character vector giving the valid types of Surv object.
# @return A Surv object.
validate_surv <- function(x, ok_types = c("right", "counting",
                                          "interval", "interval2")) {
  if (!inherits(x, "Surv"))
    stop2("LHS of 'formula' must be a 'Surv' object.")
  if (!attr(x, "type") %in% ok_types)
    stop2("Surv object type must be one of: ", comma(ok_types))
  x
}


# Return the response vector (time)
#
# @param model_frame The model frame.
# @param type The type of time variable to return:
#   "beg": the entry time for the row in the survival data,
#   "end": the exit  time for the row in the survival data,
#   "gap": the difference between entry and exit times,
#   "upp": if the row involved interval censoring, then the exit time
#          would have been the lower limit of the interval, and "upp"
#          is the upper limit of the interval.
# @return A numeric vector.
make_t <- function(model_frame, type = c("beg", "end", "gap", "upp")) {
    type <- match.arg(type)
    resp <- if (survival::is.Surv(model_frame))
                model_frame else model.response(model_frame)
    surv <- attr(resp, "type")
    err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")

    t_beg <- switch(surv,
                    "right"     = rep(0, nrow(model_frame)),
                    stop(err))

    t_end <- switch(surv,
                    "right"     = as.vector(resp[, "time"]),
                    stop(err))

    t_upp <- switch(surv,
                    "right"     = rep(NaN, nrow(model_frame)),
                    stop(err))

    switch(type,
           "beg" = t_beg,
           "end" = t_end,
           "gap" = t_end - t_beg,
           "upp" = t_upp,
           stop("Bug found: cannot handle specified 'type'."))
}

# Return the response vector (status indicator)
#
# @param model_frame The model frame.
# @return A numeric vector.
make_d <- function(model_frame) {

  resp <- if (survival::is.Surv(model_frame))
    model_frame else model.response(model_frame)
  surv <- attr(resp, "type")
  err  <- paste0("Bug found: cannot handle '", surv, "' Surv objects.")

  switch(surv,
         "right"     = as.vector(resp[, "status"]),
         stop(err))
}

parse_external <- function(external, formula, ncovs, xlevs=NULL, xbar,
                           cure_formula=NULL, ncurecovs=NULL, xlevs_cure=NULL, xbar_cure=NULL){
    ## TODO validate here.
    if (is.null(external))
        external <- list(nextern=0, stop=numeric(), start=numeric(),
                         r=integer(), n=integer(), tmax=0)
    else {
        external <- c(as.list(external), nextern=nrow(external),
                      tmax=max(external$stop))
        if (ncovs>0){
            form <- delete.response(terms(formula))
            x <- model.matrix(form, external, xlev = xlevs)
            x <- drop_intercept(x)
            x_centered <- sweep(x, 2, xbar, FUN = "-")
        }
        if (ncurecovs>0){
            form <- delete.response(terms(cure_formula))
            xcure <- model.matrix(cure_formula, external, xlev = xlevs_cure)
            xcure <- drop_intercept(xcure)
            xcure_centered <- sweep(xcure, 2, xbar_cure, FUN = "-")
        }
    }
    if ((external$nextern==0) || (ncovs==0))
        x_centered <- array(dim=c(external$nextern, ncovs))
    if ((external$nextern==0) || (ncurecovs==0))
        xcure_centered <- array(dim=c(external$nextern, ncurecovs))
    external <- c(external, nlist(x_centered), nlist(xcure_centered))
    external
}

make_basis <- function(times, basehaz, integrate = FALSE) {
  N <- length(times)
  K <- basehaz$nvars
  if (!N) { # times is NULL or empty vector
    return(matrix(0, 0, K))
  }
  mspline_basis(times, iknots = basehaz$iknots, bknots=basehaz$bknots,
                degree = basehaz$degree, integrate = integrate)
}

aa <- function(x, ...) as.array(x,...)


make_basehaz <- function(basehaz_ops,
                                times,
                                times_ext,
                                status,
                                tmin, tmax){
    df     <- basehaz_ops$df
    iknots <- basehaz_ops$iknots
    bknots <- basehaz_ops$bknots # TODO validate
    degree <- basehaz_ops$degree
    if (is.null(df))
        df <- 10L
    if (is.null(degree))
        degree <- 3L # cubic splines
    tt <- times[status == 1] # uncensored event times
    if (is.null(iknots) && !length(tt)) {
        warning2("No observed events found in the data. Censoring times will ",
                 "be used to evaluate default knot locations for splines.")
        tt <- times
    }
    ## TODO validate iknots and bknots here.  Want to allow them outside data
    if (!is.null(iknots)) {
    }
    if (is.null(bknots)){
        bknots <- c(0, tmax)
    }
    names(bknots) <- c("lower","upper")
    ttk <- unique(c(tt, times_ext))
    iknots <- get_iknots(ttk, df = df, iknots = iknots, degree = degree, intercept = TRUE)
    if (any(times<0)) stop("Some survival times are negative") # TODO should do this validation somewhere else

    nvars  <- df
    knots <- c(bknots[1], iknots, bknots[2])
    nlist(nvars, iknots, bknots, degree, df, knots)
}

# Return a vector with internal knots for 'x', based on evenly spaced quantiles
#
# @param x A numeric vector.
# @param df The degrees of freedom. If specified, then 'df - degree - intercept'.
#   knots are placed at evenly spaced percentiles of 'x'. If 'iknots' is
#   specified then 'df' is ignored.
# @param degree Non-negative integer. The degree for the spline basis.
# @param iknots Optional vector of internal knots.
# @return A numeric vector of internal knot locations, or NULL if there are
#   no internal knots.
get_iknots <- function(x, df = 5L, degree = 3L, iknots = NULL, intercept = FALSE) {

  # obtain number of internal knots
  if (is.null(iknots)) {
    nk <- df - degree - intercept
  } else {
    nk <- length(iknots)
  }
  # validate number of internal knots
    if (nk < 0) {
    stop("Number of internal knots cannot be negative.")
  }
  # if no internal knots then return empty vector
  if (nk == 0) {
    return(numeric(0))
  }
  # obtain default knot locations if necessary
  if (is.null(iknots)) {
    iknots <- qtile(x, nq = nk + 1) # evenly spaced percentiles
  }
  # return internal knot locations, ensuring they are positive
  validate_positive_scalar(iknots)

  return(iknots)
}

make_xcure <- function(cure, data){
    if (isTRUE(cure)) {
        cure_formula <- ~1
    } else if (inherits(cure,"formula")) {
        cure_formula <- cure
    } else if (isFALSE(cure)) {
        cure_formula <- NULL
    } else stop("`cure` must either be TRUE, FALSE or a model formula")
    if (!is.null(cure_formula)){
        Terms <- terms(cure_formula)
        mf <- stats::model.frame(Terms, data)
        xcure <- make_x(cure_formula, mf)
        cure <- TRUE
    } else {
        xcure <- list(ncovs = 0,
                      x_centered = matrix(nrow=nrow(data), ncol=0))
    }
    c(nlist(cure, cure_formula), xcure)
}
