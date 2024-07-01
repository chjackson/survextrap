
##' Evaluate an M-spline basis matrix at the specified times.
##'
##' Evaluate an M-spline basis matrix at the specified times.
##' Extrapolation beyond the boundary knots is done by assuming that
##' each basis term is constant beyond the boundary.
##'
##' The lower boundary is fixed to zero, and each basis term is
##' assumed to be zero at times less than zero, since these models are
##' used for hazard functions in survival data.
##'
##' @param times A numeric vector of times at which to evaluate the basis.
##'
##' @param knots Spline knots
##'
##' @param degree Spline degree
##'
##' @param integrate If \code{TRUE}, then the integrated M-spline (I-spline) basis is returned.
##'
##' @param bsmooth If \code{TRUE} then the function is constrained to
##'   also have zero derivative and second derivative at the boundary,
##'   which improves smoothness at the boundary (experimental feature).
##'
##' @references The [splines2](https://cran.r-project.org/web/packages/splines2/index.html) package is used.
##'
##' @return A two-dimensional array.  Rows are the times, and columns are the basis terms.
##'
##' @export
mspline_basis <- function(times, knots, degree=3, integrate = FALSE,
                          bsmooth = TRUE){
  if (is.null(bsmooth)) bsmooth <- TRUE
  if (bsmooth){
    if (degree != 3)
      stop("spline degree must be 3 unless using bsmooth=FALSE")
    res <- mspline_basis_bsmooth(times=times, knots=knots, integrate = integrate)
  }
  else
    res <- mspline_basis_unsmooth(times=times, knots=knots,
                                  degree=degree, integrate = integrate)
  attr(res, "times") <- times
  attr(res, "bsmooth") <- bsmooth
  attr(res, "knots") <- knots
  attr(res, "degree") <- degree
  attr(res, "Boundary.knots") <- NULL
  ## note that this overwrites the knots attributes created by splines2
  ## which are the internal knots
  res
}

mspline_basis_unsmooth <- function(times, knots, degree=3, integrate = FALSE) {
  validate_knots(knots, name="knots")
  knots <- sort(knots) # just in case
  tmax <- max(knots)
  tmin <- 0
  iknots <- knots[-length(knots)]
  bknots <- c(tmin, tmax)
  # evaluate basis at knots first, to set up use of predict()
  basis0 <- splines2::mSpline(iknots, knots = iknots, Boundary.knots = bknots,
                              degree = degree, intercept = TRUE)

  if (integrate) {
    ibasis0 <- splines2::iSpline(iknots, knots = iknots, Boundary.knots = bknots,
                                 degree = degree, intercept = TRUE)
    out <- matrix(nrow=length(times), ncol=ncol(basis0))
    iind <- times <= tmax & times >= tmin
    times_int <- times[iind]
    if (length(times_int) > 0){
        out[iind] <- predict(ibasis0, times_int)
    }
    eind <- which(times > tmax)
    ## Above the upper boundary knot
    times_ext <- times[eind]
    n_ext <- length(times_ext)
    Mmax <- predict(basis0, tmax)
    Imax <- predict(ibasis0, tmax)
    for (i in seq_len(n_ext)){
        out[eind[i],] <- Imax + Mmax*(times_ext[i] - tmax)
    }
    ## Below the lower boundary knot
    eind <- which(times < tmin & times > 0)
    times_ext <- times[eind]
    n_ext <- length(times_ext)
    Mmin <- predict(basis0, tmin)
    Imin <- predict(ibasis0, tmin)
    for (i in seq_len(n_ext)){
        out[eind[i],] <- Mmin*times_int[i]
    }
  } else {
      times <- pmin(times, tmax)
      times <- pmax(times, tmin)
      out <- predict(basis0, times)
  }
  out[times<=0,] <- 0
  aa(out)
}

##' M-spline basis that is constrained to be constant and smooth at
##' the upper boundary, so that the derivative and second derivative
##' are zero.
##'
##' Only supports cubic M-splines, i.e. degree 3
##'
##' @author Derivation by Iain Timmins (https://github.com/irtimmins)
##'
##' Lower boundary constraints not supported, on the assumption that
##' users will nearly always use a lower boundary of zero, so no need
##' to model below the boundary.
##'
##' @noRd
mspline_basis_bsmooth <- function(times, knots, integrate = FALSE) {
  basis0 <- mspline_basis_unsmooth(times, knots, degree=3, integrate)
  n <- ncol(basis0)
  res <- matrix(nrow=length(times), ncol=n-2)
  for (i in 1:(n-3))
    res[,i] <- basis0[,i]
##  res[,n-2] <- n2coef*basis0[,n-2] + n1coef*basis0[,n-1] + ncoef*basis0[,n]
  res[,n-2] <- basis0[,c(n-2,n-1,n)] %*% bsmooth_coefs(knots)
  res
}

bsmooth_coefs <- function(knots){
  knots <- sort(knots) # just in case
  iknots <- knots[-length(knots)]
  bknots <- c(0, max(knots))
  b <- splines2::mSpline(iknots, knots = iknots, Boundary.knots = bknots,
                         degree = 3, intercept = TRUE)
  b_deriv <- splines2::mSpline(iknots, knots = iknots, Boundary.knots = bknots,
                               degree = 3, intercept = TRUE, derivs=1)
  b_2deriv <- splines2::mSpline(iknots, knots = iknots, Boundary.knots = bknots,
                               degree = 3, intercept = TRUE, derivs=2)
  U <- bknots[2]
  bU <- predict(b, U)
  bdU <- predict(b_deriv, U)
  b2dU <- predict(b_2deriv, U)
  n <- length(bU)
  ncoef <- 1 / bU[n]
  n1coef <- - ncoef * bdU[n] / bdU[n-1]
  n2coef <- ( - n1coef * b2dU[n-1]  -  ncoef * b2dU[n] ) /  b2dU[n-2]
  c(n2coef, n1coef, ncoef)
}
