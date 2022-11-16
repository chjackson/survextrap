##' Obtain the cumulative hazard at a particular time, given a data frame of
##' piecewise-constant background hazards.
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
##'  The background hazard is the hazard of death from causes other than the
##'  specific cause of interest. The overall hazard is defined as the sum of the
##'  background hazard and the cause-specific hazard.  The cause-specific hazard
##'  is modelled with the M-spline model, and the background hazard is assumed
##'  to be known.
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
