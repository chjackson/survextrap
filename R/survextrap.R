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
#' @param prior_weights Spline basis weights defining the prior mean for the hazard function. By
#' default these are set to values that define a constant hazard function.  (TODO: a tool for people
#' to convert other hazard forms into prior weights). ?? prior mean or beta?
#'
#' @param cure If `TRUE` a mixture cure model is used.
#'
#' @param cure_prior Beta shape parameters for the cure probability (vector of two)
#'
#' @param basehaz_ops A list of control parameters
#'
#' `knots`: Internal knots.  If this is not supplied, then `df` chosen quantiles of the uncensored survival times.
#'
#' `bknots`: Boundary knots.  Zero and maximum time
#' ?? external?
#'
#' `df`: Degrees of freedom, i.e. the number of parameters.  Default 10.
#'
#' `degree`: polynomial degree used for the basis function (default 3, a cubic)
#'
#' @param modelid TODO
#'
#' @param algorithm TODO
#'
#' @param ... Additional arguments to supply to control the Stan fit, passed to the \pkg{rstan} function \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}}
#'
#'
#' @export
survextrap <- function(formula,
                       data,
                       external=NULL,
                       smooth_sd = "bayes",
                       prior_weights = NULL,
                       cure = FALSE,
                       cure_prior = c(1,1),
                       basehaz_ops = NULL,
                       modelid = 1,
                       algorithm = "mcmc",
                       ...)
{

    mftmp <- make_model_frame(formula, data)
    mf <- mftmp$mf # model frame

    t_beg <- make_t(mf, type = "beg") # entry time
    t_end <- make_t(mf, type = "end") # exit  time
    t_upp <- make_t(mf, type = "upp") # upper time for interval censoring
    status <- make_d(mf)
    t_event <- aa(t_end[status == 1]) # exact event time
    t_rcens <- aa(t_end[status == 0]) # right censoring time
    nevent <- sum(status == 1)
    nrcens <- sum(status == 0)
    ind_event <- which(status==1)
    ind_rcens <- which(status==0)

    make_x <- function(formula, mf, xlevs=NULL){
        x <- model.matrix(formula, mf, xlev = xlevs)
        x <- drop_intercept(x)
        ncovs <- NCOL(x)
        xbar <- aa(colMeans(x))
        x_centered <- sweep(x, 2, xbar, FUN = "-")

        mf2 <- mf[-attr(terms(mf), "response")]
        if (ncol(mf2) > 0){
            xinds <- list(
                factor = which(sapply(mf2, is.factor)),
                numeric = which(sapply(mf2, is.numeric))
            )
            xlevs <- lapply(mf2[xinds$factor], levels)
            if ((length(xinds$factor)==1) && (length(xinds$numeric)==0)){
                base_levs <- xlevs
            }
            else if (length(xinds$factor) > 0)
                base_levs <- lapply(xlevs, function(x)x[1])
            else base_levs <- NULL
            mean_of_numerics <- lapply(mf2[xinds$numeric], mean)
            mf_baseline <- as.data.frame(c(mean_of_numerics, base_levs))
        } else xlevs <- mf_baseline <- xinds <- NULL

        nlist(x, x_centered, xbar, N = NROW(x),
              ncovs = NCOL(x),
              xnames=colnames(x),
              xlevs, xinds,
              mf_baseline)
    }

    x <- make_x(formula, mf)
    ncovs <- x$ncovs
    x_event <- x$x_centered[ind_event, , drop = FALSE]
    x_rcens <- x$x_centered[ind_rcens, , drop = FALSE]

    t_tmp <- sum(rowMeans(cbind(t_end, t_upp), na.rm = TRUE) - t_beg)
    d_tmp <- sum(!status == 0)
    log_crude_event_rate <- log(d_tmp / t_tmp)
    if (is.infinite(log_crude_event_rate))
        log_crude_event_rate <- 0 # avoids error when there are zero events

    external <- parse_external(external, formula, ncovs, xlevs=x$xlevs, xbar=x$xbar)
    t_ext_stop <- aa(external$stop)
    t_ext_start <- aa(external$start)
    r_ext <- aa(external$r)
    n_ext <- aa(external$n)
    nextern <- external$nextern
    x_ext <- external$x_centered

    basehaz <- handle_basehaz_surv(basehaz_ops    = basehaz_ops,
                                   times          = t_end,
                                   times_ext      = unique(c(t_ext_start, t_ext_stop)),
                                   status         = status,
                                   min_t          = min(t_beg),
                                   max_t          = max(c(t_end,t_upp,external$tmax), na.rm = TRUE))

    basis_event  <- make_basis(t_event, basehaz)
    ibasis_event <- make_basis(t_event, basehaz, integrate = TRUE)
    ibasis_rcens <- make_basis(t_rcens, basehaz, integrate = TRUE)
    nvars <- basehaz$nvars

    if (is.null(prior_weights)){
        prior_weights <- mspline_uniform_weights(knots = basehaz$iknots, bknots=basehaz$bknots)
    } else {
        ## TODO validate prior weights, reconcile with betaraw
    }
    beta_mean <- log(prior_weights[-1] / prior_weights[1])
    cure_prob_init <- if (cure) 0.5 else aa(numeric())
    if (!(modelid %in% 1:2)) stop("modelid should be 1 or 2")

    if (modelid==2){
        ## nvars is number of parameters including the intercept
        nvars <- 2
        ## "basis" is the observed times and "ibasis" is the censored
        basis_event <- array(t_event, dim=c(length(t_event), 2)) # just first col of these used
        ibasis_event <- array(t_event, dim=c(length(t_event), 2))
        ibasis_rcens <- array(t_rcens, dim=c(length(t_rcens), 2))
        beta_mean <- aa(1) # prior mean for log shape parameter. shape is coefs[1], log scale is eta[1]
    }

    ibasis_ext_stop <- if (nextern>0) make_basis(t_ext_stop, basehaz, integrate = TRUE) else matrix(nrow=0, ncol=nvars)
    ibasis_ext_start <- if (nextern>0) make_basis(t_ext_start, basehaz, integrate = TRUE) else matrix(nrow=0, ncol=nvars)

    est_smooth <- (smooth_sd == "bayes")
    if (est_smooth) smooth_sd <- 1
    standata <- nlist(nevent, nrcens, nvars, nextern, ncovs,
                      log_crude_event_rate,
                      basis_event, ibasis_event, ibasis_rcens,
                      ibasis_ext_stop, ibasis_ext_start,
                      x_event, x_rcens,
                      r_ext, n_ext, x_ext,
                      beta_mean,
                      est_smooth,
                      cure,
                      cure_shape = cure_prior, ## TODO standardise name, mean+ESS param
                      modelid
                      )
    staninit <- list(gamma = aa(0),
                     loghr = aa(rep(0, ncovs)),
                     beta_err = rep(0, nvars),
                     smooth_sd = aa(smooth_sd),
                     cure_prob = cure_prob_init)
    if (identical(smooth_sd, "eb")){
        smooth_sd <- eb_smoothness(standata, staninit)
    }
    standata$smooth_sd_fixed <- if (est_smooth) aa(numeric()) else aa(smooth_sd)

    if (algorithm=="opt")
        fits <- rstan::optimizing(stanmodels$survextrap, data=standata, init=staninit)
    else 
        fits <- rstan::sampling(stanmodels$survextrap, data=standata, 
                                pars = "beta", include=FALSE, 
                                ...)

    km <- survminer::surv_summary(survfit(formula, data=data), data=data)

    res <- list(formula=formula,
                stanfit=fits, standata=standata, basehaz=basehaz,
                entrytime=t_beg, eventtime=t_end, external=external, # TODO way to strip the data
                nvars = nvars,
                ncovs = ncovs,
                prior_weights = prior_weights,
                smooth_sd = smooth_sd,
                cure = cure,
                modelid = modelid,
                km = km,
                xnames = x$xnames,
                xlevs = x$xlevs,
                xinds = x$xinds,
                xbar = x$xbar,
                log_crude_event_rate = log_crude_event_rate,
                mf_baseline = x$mf_baseline)
    class(res) <- "survextrap"
    res
}


eb_smoothness <- function(standata, staninit){
    standata$est_smooth <- 1
    standata$smooth_sd_fixed  <- aa(1) # dummy
    fits <- rstan::optimizing(stanmodels$survextrap, data=standata, init=staninit, hessian=FALSE, verbose=TRUE)
    if (fits$return_code==0){
        smooth_sd <- fits$par["smooth_sd[1]"]
    } else {
        warning("Empirical Bayes estimation of smoothness parameter failed, continuing with default")
        smooth_sd <- 1
    }
    smooth_sd
}

parse_formula_and_data <- function(formula, data) {
  formula <- validate_formula(formula, needs_response = TRUE)
  allvars <- all.vars(formula)
  allvars_form <- reformulate(allvars)

  # LHS of entire formula
  lhs       <- lhs(formula)         # LHS as expression
  lhs_form  <- reformulate_lhs(lhs) # LHS as formula

  # RHS of entire formula
  rhs       <- rhs(formula)         # RHS as expression
  rhs_form  <- reformulate_rhs(rhs) # RHS as formula

  # evaluate model data (row subsetting etc)
  data <- make_model_data(allvars_form, data)

  # evaluated response variables
  surv <- eval(lhs, envir = data) # Surv object
  surv <- validate_surv(surv)
  type <- attr(surv, "type")

  if (type == "right") {
    min_t    <- 0
    max_t    <- max(surv[, "time"])
    status   <- as.vector(surv[, "status"])
    t_end    <- as.vector(surv[, "time"])
  }  else {
      stop("Right-censoring is only kind of censoring supported")
  }

  if (any(is.na(status)))
    stop2("Invalid status indicator in Surv object.")

  if (any(status < 0 || status > 3))
    stop2("Invalid status indicator in Surv object.")

  type <- attr(surv, "type")

  nlist(formula,
        data,
        allvars,
        allvars_form,
        lhs,
        lhs_form,
        rhs,
        rhs_form,
        surv_type = attr(surv, "type"))

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


## omitted xlevs, drop.unused.levels

make_model_frame <- function(formula,
                             data,
                             na.action){
  Terms <- terms(formula) # no random effects for now
  mf <- stats::model.frame(Terms,
                           data,
                           na.action = na.action)
    list(mf=mf,
         mt=attr(mf, "Terms"))
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

parse_external <- function(external, formula, ncovs, xlevs=NULL, xbar){
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
    }
    if ((external$nextern==0) || (ncovs==0))
        x_centered <- array(dim=c(external$nextern, ncovs))
    external <- c(external, nlist(x_centered))
    external
}

make_basis <- function(times, basehaz, integrate = FALSE) {
  N <- length(times)
  K <- basehaz$nvars
  if (!N) { # times is NULL or empty vector
    return(matrix(0, 0, K))
  }
  switch(basehaz$type_name,
         "ms"          = mspline_basis(times, iknots = basehaz$iknots, bknots=basehaz$bknots,
                                       degree = basehaz$degree, integrate = integrate),
         stop2("Bug found: unknown type of baseline hazard."))
}

aa <- function(x, ...) as.array(x,...)


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
    ## FIXME why does this not print the error sometines? 
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


handle_basehaz_surv <- function(basehaz_ops,
                                times,
                                times_ext,
                                status,
                                min_t, max_t){

    df     <- basehaz_ops$df
    knots  <- basehaz_ops$knots
    degree <- basehaz_ops$degree
    if (is.null(df))
        df <- 10L
    if (is.null(degree))
        degree <- 3L # cubic splines
    tt <- times[status == 1] # uncensored event times
    if (is.null(knots) && !length(tt)) {
        warning2("No observed events found in the data. Censoring times will ",
                 "be used to evaluate default knot locations for splines.")
        tt <- times
    }
    ## TODO validate knots here.  Want to allow them outside data
    if (!is.null(knots)) {
    }
    bknots <- basehaz_ops$bknots # TODO validate
    if (is.null(bknots)){
        bknots <- c(0, max_t)
    }
    ttk <- unique(c(tt, times_ext))
    iknots <- get_iknots(ttk, df = df, iknots = knots, degree = degree, intercept = TRUE)
    if (any(times<0)) stop("Some survival times are negative") # TODO should do this validation somewhere else

    nvars  <- df

    nlist(type_name = "ms",
          nvars,
          iknots,
          bknots,
          degree,
          df = nvars,
          user_df = nvars,
          knots = c(bknots[1], iknots, bknots[2]))

}

