
##' Create a default M-spline model structure
##'
##' @param df Desired number of basis terms, or "degrees of freedom"
##'   in the spline.  If \code{knots} is not supplied, the number of
##'   knots is then chosen to satisfy this.
##' 
##' @param degree Spline polynomial degree.  Can only be changed from
##' the default of 3 if \code{bsmooth} is \code{FALSE}.
##' 
##' @param bsmooth If \code{TRUE} then the function is constrained to
##'   also have zero derivative and second derivative at the boundary.
##'
##' @param knots Vector of knot locations. If not supplied, \code{df}
#'   has to be specified.  One of two rules is then used to choose the
#'   knot locations.  If \code{bknots} is specified, a set of equally
#'   spaced knots between zero and \code{bknots} is used.  Otherwise
#'   if \code{obstimes} is supplied, the knots are chosen as equally
#'   spaced quantiles of \code{obstimes}.
#'
#'   The number of knots (excluding zero) is \code{df - degree + 1} if
#'   \code{bsmooth} is \code{TRUE}, or \code{df - degree - 1}
#'   otherwise.
##' 
##' @param bknot Location of the final spline knot.
##' 
##' @param obstimes Vector of observation times whose quantiles will be
##' used to choose knot locations
##'
##' @return A list with fundamental components \code{knots},
##'   \code{degree}, and \code{bsmooth}, as documented above.
##'
##' The component \code{df} is also included, and derived as a consequence
##' of the fundamental components.
##'
##' @export
mspline_init <- function(df = 10,
                         degree = 3,
                         bsmooth = TRUE,
                         knots = NULL,
                         bknot = 10,
                         obstimes = NULL
                         )
{
  if (!is.numeric(bknot) || (bknot <=0))
    stop("`bknot` should be a positive number")
  
  if (is.null(bsmooth))
    bsmooth <- TRUE
  else if (!is.logical(bsmooth))
    stop("bsmooth should be NULL, TRUE or FALSE")
  if (is.null(degree))
    degree <- 3
  else if (degree!=3 && bsmooth)
    stop("M-spline degree must be 3 unless bsmooth=FALSE")

  knots <- sort(knots)
  if (!is.null(knots)) {
    validate_knots(knots, "knots")
    if (bsmooth) 
      df <- length(knots) + degree - 2
    else df <- length(knots) + degree
  }
  if (is.null(df))
    df <- 10
  else 
    validate_df(df, degree, bsmooth)

  if (is.null(knots)){
    nk <- if (bsmooth) df - degree + 2 else df - degree
    if (!is.null(obstimes)){
      knots <- quantile(obstimes, probs=seq(1, nk)/nk)
    } else if (!is.null(bknot)) {
      knots <- seq(0, bknot, length.out=nk+1)[-1]
    }
    else stop("`knots` not supplied, and no way to choose default knots")
  }
  
  nlist(knots, degree, bsmooth, df)
}

##' Validate an M-spline object supplied as a list, choosing defaults
##' if needed.
##'
##' @param mspline A list with any or none of the following components:
##' \code{df}, \code{degree}, \code{bsmooth}, \code{knots}, \code{bknot},
##' as documented in \code{\link{mspline_init}}.
##'
##' @inheritParams mspline_init
##'
##' @return A list defining the M-spline, with any list components set
##'   to defaults.  See \code{\link{mspline_init}} for details.
##'
##' If \code{mspline$knots} is not supplied, giving knot locations, then
##' either \code{mspline$bknot} or \code{obstimes} must be specified,
##' so that default locations can be obtained. 
##'
##' @export
mspline_list_init <- function(mspline, obstimes=NULL){
  if (!is.null(mspline)){
    if (!is.list(mspline)) stop("`mspline` should be a list")
    bh_names <- c("df","knots","degree","nvars","knots","bsmooth","add_knots")
    bad_names <- setdiff(names(mspline), bh_names)
    if (length(bad_names) > 0) {
      blist <- paste(bad_names, collapse=",")
      plural <- if (length(bad_names) > 1) "s" else ""
      warning(sprintf("Element%s `%s` of `mspline` is unused. Ignoring.", plural, blist))
    }
  } else mspline <- list()
  mspline <- mspline_init(df = mspline$df, degree = mspline$degree,
                          bsmooth = mspline$bsmooth, knots = mspline$knots,
                          obstimes = obstimes)
  mspline
}

##' Create an M-spline survival model, both structure and parameters.
##'
##' \code{\link{mspline_init}} is first used to create the M-spline
##' model structure, including knot positions.  Parameters including
##' basis coefficients and scale are either supplied or set to a
##' default that defines a constant hazard model.
##'
##' This function is not for fitting models to data, but for setting
##' up a theoretical M-spline model for illustration.
##'
##' @inheritParams mspline_init
##' 
##' @param coefs Basis coefficients 
##'
##' @param hscale Hazard scale parameter 
##'
##' @export
msplinemodel_init <- function(df = 10,
                              degree = 3,
                              bsmooth = TRUE,
                              knots = NULL,
                              bknot = 10,
                              obstimes = NULL,
                              coefs = NULL,
                              hscale = 1)
{
  mspline <- mspline_init(df=df, degree=degree,
                          bsmooth=bsmooth, knots=knots,
                          bknot=bknot, obstimes=obstimes)
  if (is.null(coefs))
    mspline$coefs <- mspline_constant_coefs(mspline)
  else 
    validate_coefs(coefs, nvars=mspline$df)
  validate_hscale(hscale)
  mspline$hscale <- hscale
  mspline
}
