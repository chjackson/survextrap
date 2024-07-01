
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
#'   knot locations.  If \code{bknot} is specified, a set of equally
#'   spaced knots between zero and \code{bknot} is used.  Otherwise
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
  else if (degree < 0 || !is.numeric(degree))
    stop("M-spline degree must be a nonnegative number")

  knots <- sort(knots)
  if (!is.null(knots)) {
    validate_knots(knots, "knots")
    df <- mspline_df(knots, degree, bsmooth)
  }
  if (is.null(df))
    df <- 10
  else
    validate_df(df, degree, bsmooth)

  if (is.null(knots)){
    nk <- if (bsmooth) df - degree + 2 else df - degree
    if (!is.null(obstimes)){
      knots <- quantile(obstimes, probs=seq(1, nk)/nk)
      names(knots) <- round(knots, 1)
    } else if (!is.null(bknot)) {
      knots <- seq(0, bknot, length.out=nk+1)[-1]
    }
    else stop("`knots` not supplied, and no way to choose default knots")
  }

  mspline_update(nlist(knots, degree, bsmooth))
}

mspline_df <- function(knots, degree, bsmooth){
  length(knots) + degree - 2*bsmooth
}

#' Deduce characteristics of an M-spline that are defined deterministically
#' given the knots and degree (and whether we want to smooth the upper boundary)
#'
#' @noRd
mspline_update <- function(mspline){
  if (is.null(mspline$bsmooth)) mspline$bsmooth <- TRUE
  mspline$df <- mspline_df(mspline$knots, mspline$degree, mspline$bsmooth)
  mspline$basis_means <- make_basis_means(mspline$knots, mspline$degree, mspline$bsmooth)
  mspline$basis_spans <- make_basis_spans(mspline$knots, mspline$degree, mspline$bsmooth)
  mspline$sqrt_wt <- sqrt(mspline$basis_spans / sum(mspline$basis_spans))
  mspline
}

#' Mean of a random variable whose distribution is defined by a single
#' M-spline basis term.  This mean describes where the basis term is
#' "centred".  Used, e.g. to construct the random walk prior for basis
#' coefficients.
#'
#' @noRd
make_basis_means <- function(knots, degree, bsmooth){
  df <- mspline_df(knots, degree, bsmooth)
  order <- degree + 1
  lower_boundary <- 0
  knots_expand <- c(rep(lower_boundary, order), knots, rep(max(knots), order-1))
  basis_means <- numeric(df)
  df_unsmooth <- length(knots) + degree
  bmtmp <- numeric(df_unsmooth)
  for (i in 1:df_unsmooth){
    bmtmp[i] <- sum(knots_expand[i:(i+order)]) / (order + 1)
  }
  if (bsmooth){
    for (i in 1:(df-1))
      basis_means[i] <- bmtmp[i]
    basis_means[df] <- bmtmp[df_unsmooth - c(2,1,0)] %*% bsmooth_coefs(knots)
    if (basis_means[df] < basis_means[df-1])
      basis_means[df] <- bmtmp[df] # arbitrarily
  } else {
    basis_means <- bmtmp
  }
  basis_means
}

#' Determine the differences between lag-3 knots from the second to
#' the second last knot.  Used to construct the weights for the random
#' walk prior.
#'
#' From Li and Cao https://arxiv.org/abs/2201.06808
#' 
#' @noRd
make_basis_spans <- function(knots, degree, bsmooth){
  df <- mspline_df(knots, degree, bsmooth)
  order <- degree + 1
  lower_boundary <- 0
  knots_expand <- c(rep(lower_boundary, order), knots, rep(max(knots), order-1))
  if (order==1)
    basis_spans <- knots_expand[2:df] - knots_expand[1:(df-1)]
  else 
    basis_spans <- knots_expand[(order+1):(df+order-1)] - knots_expand[2:df]
  basis_spans
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
##' @return A list defining the M-spline, with any omitted list components set
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
    bh_names <- c("df","knots","degree","nvars","knots","bsmooth","add_knots","basis_means","basis_spans","sqrt_wt")
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
##' @return A list defining the M-spline, with any omitted list
##'   components set to defaults.  See \code{\link{mspline_init}} for
##'   details.  The parameters are included as the \code{coefs} and
##'   \code{hscale} components.
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

##' Determine M-spline basis coefficients which give a constant function.
##'
##' This works by obtaining the coefficients of the corresponding
##' B-spline basis, which are equal if the B-spline is a constant
##' function.
##'
##' @param mspline A list with components `knots` (vector of knots),
##' `degree` (polynomial degree) and `bsmooth` (logical for smoothness
##' constraint at boundary), defining an M-spline configuration.
##'
##' @param logit If \code{TRUE} then the multinomial logit transform of the coefficients
##' is returned.  This is a vector of length one less than the number of coefficients,
##' with the rth element defined by \eqn{log(coefs[r+1] / coefs[1])}.
##'
##' @references Ramsay, J. O. (1988). Monotone regression splines in action. Statistical Science, 3(4): 425-441.
##'
##' @export
mspline_constant_coefs <- function(mspline, logit=FALSE){
  mspline <- mspline_list_init(mspline)
  iknots <- mspline$knots[-length(mspline$knots)]
  bknots <- c(0, max(mspline$knots))
  degree <- mspline$degree

  ## Firstly determine coefs for constant function under unsmoothed basis
  knot_seq <- c(rep(bknots[1], degree+1), iknots, rep(bknots[2], degree+1))
  K <- length(iknots) + degree + 1
  p_const <- numeric(K)
  rescale <- (degree+1)*(bknots[2] - bknots[1])
  for(i in 1:K){
    p_const[i] <- (knot_seq[i+degree+1] - knot_seq[i]) / rescale
  }

  ## Deduce equivalent coefs for bsmoothed basis with same knots
  if (mspline$bsmooth)
    p_const <- c(p_const[1:(K-3)],
                 p_const[K]*mspline_basis_unsmooth(bknots[2], mspline$knots)[K])

  p_const <- p_const/sum(p_const)

  if (logit) log(p_const[-1]/p_const[1]) else p_const
}
