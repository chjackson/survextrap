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
##' @param strata A vector of the same length as \code{t} indicating strata.
##' There should then be
##'
##' * a column in \code{backhaz} named \code{stratum}, which should
##' be a factor with the same levels as \code{strata} (or analogous
##' character vector).
##'
##' * a row in \code{backhaz} indicating the background hazard at each
##' time interval for each stratum in \code{strata}.
##'
##' We assume that the background hazard is known at all times, and defined by
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
##' @md
##' @keywords internal
get_cum_backhaz <- function(t, backhaz, strata=NULL){
  res <- numeric(length(t))
  if (is.null(strata)) {
    strata <- rep(1, length(t))
    backhaz$stratum <- 1
  } else {
    if (is.null(backhaz$stratum))
      stop("`strata` provided, but no column named `stratum` in `backhaz`")
  }
  nst <- length(unique(strata))
  for (j in 1:nst){
    stj <- unique(strata)[j]
    ind <- backhaz$stratum==stj
    if (sum(ind)==0) stop(sprintf("No stratum with value %s found in `backhaz`",
                                  stj))
    backhazj <- backhaz[ind,,drop=FALSE]
    validate_backhaz_onestratum(backhazj, stj, nst)
    time <- backhazj$time
    haz <- backhazj$hazard
    timelag <- diff(time)
    tj <- t[strata==stj]
    ind <- findInterval(tj, time)
    cumhaz <- haz[ind]*(tj - time[ind])
    for (i in seq_along(tj)){
      if (ind[i] > 1)
        cumhaz[i] <- cumhaz[i] + sum(timelag[1:(ind[i]-1)] * haz[1:(ind[i]-1)])
    }
    res[strata==stj] <- cumhaz
  }
  res
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
  if (any(backhaz$hazard < 0))
    stop("backhaz$hazard should all be non-negative")
}

validate_backhaz_onestratum <- function(backhaz, stratum, nst){
  str <- if (nst==1) "" else sprintf(" within stratum \"%s\"", stratum)
  if (backhaz$time[1] != 0)
    stop(sprintf("The first element of the \"time\" column of `backhaz` %s should be 0", str))
  if (is.unsorted(backhaz$time))
    stop(sprintf("backhaz$time should be sorted in increasing order%s", str))
}

#' Parse the background hazard data supplied to survextrap, and return
#' the essential information needed for the likelihood computation.
#'
#' @inheritParams survextrap
#'
#' @param td Output of \code{make_td} giving observed event and censoring
#' times in the individual data
#'
#' @return List with components:
#' 
#' \code{event} Background hazard at observed event times.  This was
#' either supplied to \code{survextrap} directly, or via a data frame
#' of piecewise constant hazards and their change times, potentially
#' stratified by factors in the individual data.
#'
#' \code{df} Data frame giving background hazard at all times (if
#' \code{backhaz} was supplied to \code{survextrap} in this form)
#'
#' \code{strata} Background hazard strata, as supplied in
#' the \code{backhaz_strata} argument to \code{survextrap}.
#'
#' \code{relative} Indicator for whether a relative survival (i.e. an
#' additive hazards) model is being used, i.e. whether \code{backhaz} was
#' supplied to \code{survextrap}
#'
#' @noRd 
make_backhaz <- function(backhaz, data, external, td, backhaz_strata){
  if (is.data.frame(backhaz)) {
    if (!is.null(data)){
      st <- make_backhaz_strata(backhaz_strata, backhaz, data)
      backhaz_event <- numeric(td$nevent)
      for (i in 1:st$nst){
        bhi <- backhaz[st$back==st$strata[i],,drop=FALSE]
        validate_backhaz_onestratum(bhi, st$strata[i], st$nst)
        eventtimesi <- td$t_event[st$dat[td$ind_event] == st$strata[i]]
        backhaz_event[st$dat[td$ind_event]==st$strata[i]] <-
          bhi$hazard[findInterval(eventtimesi, bhi$time)]
      }
      backhaz_df <- backhaz
      backhaz_df$stratum <- factor(st$back)
    } else if (is.data.frame(external)) { ## if no IPD
      st <- make_backhaz_strata(backhaz_strata, backhaz, external)
      backhaz_df <- backhaz
      backhaz_event <- numeric(0)
      backhaz_df$stratum <- factor(st$back)
    } else { # external data in invalid form, not yet validated at this point
      backhaz_event <- numeric(0)
      backhaz_df <- backhaz
    }
  } else if (!is.null(backhaz)) {
    if (!is.character(backhaz))
      stop("`backhaz` should be a data frame or a character string giving the name of a column in the data")
    if (!is.null(backhaz_strata))
      warning("Ignoring `backhaz_strata`, which requires `backhaz` to be supplied as a data frame")
    bh <- data[[backhaz]]
    backhaz_event <- bh[td$ind_event]
    backhaz_df <- NULL
  } else backhaz_event <- backhaz_df <- NULL
  relative <- !is.null(backhaz)
  backhaz_event <- if (relative) aa(backhaz_event) else aa(numeric(td$nevent))
  list(event = backhaz_event, df=backhaz_df,
       strata=backhaz_strata, relative=relative)
}



#' Check strata specification is valid: that indicated strata are present in the data.
#'
#' @param strata Character vector identifying strata, or \code{NULL} if no strata.  These should
#' correspond to names of variables in both \code{backhaz} and \code{data}.
#'
#' @param data Original individual-level dataset
#'
#' @param backhaz Background hazard data frame.  Not relevant when backhaz is a variable in \code{data},
#' as we don't need to match it by strata in this case.
#'
#' @param external Original external dataset
#'
#' @noRd
validate_backhaz_strata <- function(strata=NULL, data, backhaz, external=NULL){
  if (is.null(strata)) return()
  if (!is.character(strata))
    stop("`strata` should be a character vector of variables indicating background hazard strata")
  check_strata_in(strata, backhaz)
  check_strata_in(strata, data)
  check_strata_in(strata, external)
}

check_strata_in <- function(strata, dat){
  if (is.null(dat)) return()
  bad_strata <- strata[!(strata %in% names(dat))]
  plural <- if (length(bad_strata) > 1) "s" else ""
  if (length(bad_strata) > 0)
    stop(sprintf("strata variable%s %s not found in `%s`",
                 plural, paste(sprintf("\"%s\"",bad_strata),collapse=","),
                 deparse(substitute(dat))))
}

#' Standardise strata variables between data and background hazards,
#' by checking they have a consistent type between the datasets, and
#' all data can be matched to a background hazard.  Return a single
#' stratum variable formed by pasting together all stratifying factors.
#'
#' @param dat either the individual or external data
#'
#' @return A list with components:
#'
#' \code{dat} vector with with one element for each row of \code{dat},
#' giving the strata.
#'
#' \code{back} vector with with one element for each row of
#' \code{backhaz}, giving the strata.
#'
#' \code{strata} unique values of the strata in \code{dat}.
#'
#' \code{nst} number of unique strata.
#'
#' If \code{backhaz_strata=NULL}, one single stratum with value 1
#' is created, and the list is filled in accordingly.
#'
#' @noRd
make_backhaz_strata <- function(backhaz_strata, backhaz, dat){
  if (is.null(backhaz)) return(NULL)
  st <- backhaz_strata
  if (is.null(dat))
    return(list(
      dat = 1, back = rep(1, nrow(backhaz)),
      strata = 1,  nst=1
    ))
  else if (is.null(st))
    return(list(
      dat = rep(1, nrow(dat)),
      back = rep(1, nrow(backhaz)),
      strata = 1, nst=1
    ))
  validate_backhaz_df(backhaz)
  validate_backhaz_strata(backhaz_strata, dat, backhaz)
  for (i in seq_along(st)){
    ## check mismatch between factor/character and integer/numeric
    dclass <- class(dat[[st[i]]])
    bclass <- class(backhaz[[st[i]]])
    ## can match character to factor, and integer to numeric
    ## help(match): factors, raw vectors and lists are converted to character vectors
    if (((dclass %in% c("factor", "character")) && (!(bclass %in% c("factor", "character")))) ||
        ((bclass %in% c("factor", "character")) && (!(dclass %in% c("factor", "character")))))
      stop(sprintf("stratifying variable \"%s\" has class \"%s\" in `%s`, but class \"%s\" in `backhaz`, so cannot match",
                   st[i], dclass, deparse(substitute(dat)), bclass))
    firstmatchi <- match(dat[[st[i]]], backhaz[[st[i]]])
    badindiv <- which(is.na(firstmatchi))
    if (length(badindiv) > 0){
      badindiv_str <- if (length(badindiv) > 3)
                        paste(paste(badindiv[1:3],collapse=","), "and others")
                      else paste(badindiv,collapse=",")
      stop(sprintf("rows %s in `%s` have no `%s` value in `backhaz` matching their value for `%s` in `%s`",
                   badindiv_str, deparse(substitute(dat)),
                   st[i], st[i], deparse(substitute(dat))))
    }
  }
  strata_dat <- do.call(function(...)paste(...,sep=";"), dat[,st,drop=FALSE])
  strata_back <- do.call(function(...)paste(...,sep=";"), backhaz[,st,drop=FALSE])
  list(dat = strata_dat,
       back = strata_back,
       strata = unique(strata_dat),
       nst = length(unique(strata_dat)))
}
