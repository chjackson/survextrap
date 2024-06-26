% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mspline_plots.R
\name{mspline_plotsetup}
\alias{mspline_plotsetup}
\title{Get basis for an illustration of an M-spline with given knots.}
\usage{
mspline_plotsetup(
  knots,
  bknot = 10,
  tmin = NULL,
  tmax = NULL,
  degree = 3,
  df = 10,
  bsmooth = TRUE
)
}
\arguments{
\item{knots}{Vector of knot locations. If not supplied, \code{df}
has to be specified.  One of two rules is then used to choose the
knot locations.  If \code{bknot} is specified, a set of equally
spaced knots between zero and \code{bknot} is used.  Otherwise
if \code{obstimes} is supplied, the knots are chosen as equally
spaced quantiles of \code{obstimes}.

The number of knots (excluding zero) is \code{df - degree + 1} if
\code{bsmooth} is \code{TRUE}, or \code{df - degree - 1}
otherwise.}

\item{bknot}{Location of the final spline knot.}

\item{tmin}{Minimum plotting time.  Defaults to zero.}

\item{tmax}{Maximum plotting time.  Defaults to the highest knot.}

\item{degree}{Spline polynomial degree.  Can only be changed from
the default of 3 if \code{bsmooth} is \code{FALSE}.}

\item{df}{Desired number of basis terms, or "degrees of freedom"
in the spline.  If \code{knots} is not supplied, the number of
knots is then chosen to satisfy this.}

\item{bsmooth}{If \code{TRUE} then the function is constrained to
also have zero derivative and second derivative at the boundary.}
}
\value{
Data frame containing the basis, as returned by \code{\link{mspline_basis}}.
}
\description{
Get basis for an illustration of an M-spline with given knots.
}
