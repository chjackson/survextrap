##' Cumulative background hazards
##'
##' Obtain the cumulative hazard at a particular time, given a data
##' frame of piecewise-constant background hazards.  This is used in
##' the computation of likelihoods and results in additive hazards
##' models.
##'
##' @param t A time point, or vector of time points, to calculate the cumulative
##'  hazard at.
##'
##' @param backhaz A data frame with two columns, named \code{"time"} and
##'  \code{"hazard"}. Each row gives the "background" hazard between the
##'  specified time and the next time.
##'
##'  We assume that the background hazard is known at all times, and defined by
##'  a piecewise-constant function (step function) of time. The first element of
##'  \code{"time"} should be 0, and the final row specifies the hazard at all
##'  times greater than the last element of \code{"time"}.
##'
##' In additive hazards models, the background hazard is the hazard of
##'  death from causes other than the specific cause of interest. The
##'  overall hazard is defined as the sum of the background hazard and
##'  the cause-specific hazard.  In \code{\link{survextrap}}, the
##'  cause-specific hazard is modelled with the M-spline model, and
##'  the background hazard is assumed to be known.
##'
##' @return The cumulative hazard at `t`.
##'
##' @keywords internal
get_cum_backhaz <- function(t, backhaz){
  time <- backhaz$time
  haz <- backhaz$hazard
  timelag <- diff(time)
  ind <- findInterval(t, time)
  cumhaz <- haz[ind]*(t - time[ind])
  for (i in seq_along(t)){
    if (ind[i] > 1)
      cumhaz[i] <- cumhaz[i] + sum(timelag[1:(ind[i]-1)] * haz[1:(ind[i]-1)])
  }
  cumhaz
}

validate_backhaz_df <- function(backhaz){
  if (!is.data.frame(backhaz))
    stop("`backhaz` should be a data frame")
  if (!("time" %in% names(backhaz)))
    stop("`backhaz` should have a column called \"time\"")
  if (!("hazard" %in% names(backhaz)))
    stop("`backhaz` should have a column called \"hazard\"")
  if (!is.numeric(backhaz$hazard))
    stop("backhaz$hazard should be numeric")
  if (!is.numeric(backhaz$time))
    stop("backhaz$time should be numeric")
  if (backhaz$time[1] != 0)
    stop("The first element of the \"time\" column of `backhaz` should be 0")
  if (is.unsorted(backhaz$time))
    stop("backhaz$time should be sorted in increasing order")
  if (any(backhaz$hazard < 0))
    stop("backhaz$hazard should all be non-negative")
}


make_backhaz <- function(backhaz, data, td){
  if (is.data.frame(backhaz)) {
    validate_backhaz_df(backhaz)
    backhaz_event <- backhaz$hazard[findInterval(td$t_event, backhaz$time)]
    backhaz_df <- backhaz
  } else if (!is.null(backhaz)) {
    if (!is.character(backhaz)) stop("`backhaz` should be a data frame or a character string giving the name of a column in the data")
    bh <- data[[backhaz]]
    backhaz_event <- bh[td$ind_event]
    backhaz_df <- NULL
  } else backhaz_event <- backhaz_df <- NULL
  relative <- !is.null(backhaz_event)
  backhaz_event <- if (relative) aa(backhaz_event) else aa(numeric(td$nevent))
  list(event = backhaz_event, df=backhaz_df, relative=relative)
}
