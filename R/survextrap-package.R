#' The 'survextrap' package.
#'
#' @description For an introduction to and overview of the `survextrap` package, and full documentation, see
#'
#' \url{https://chjackson.github.io/survextrap/}.
#'
#' For details of the methods, see the paper by Jackson (2023). 
#' 
#' @references Jackson, C. (2023) \code{survextrap}: a package for flexible and transparent
#' survival extrapolation. BMC Medical Research Methodology 23:282.
#' \doi{10.1186/s12874-023-02094-1}
#' 
#' Timmins I, Torabi F, Jackson C, Lambert P, Sweeting M J. (2025)
#' Simulation-based assessment of a Bayesian survival model with flexible baseline hazard and time-dependent effects.
#' \doi{10.48550/arXiv.2503.21388}.
#' 
#' @name survextrap-package
#' @useDynLib survextrap, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling nlist
#' @import ggplot2
#' @importFrom survival survfit
#' @importFrom posterior as_draws as_draws_matrix rhat ess_bulk
#' @importFrom gridExtra grid.arrange
#' @importFrom splines2 mSpline iSpline
#' 
#' @importFrom stats as.formula delete.response dweibull integrate median model.matrix model.response optim predict pweibull quantile qbeta qgamma qnorm qt rlogis rnorm runif rbeta rgamma rt sd terms time var formula model.frame na.pass plogis qlogis get_all_vars reshape .getXlevels setNames
#'
"_PACKAGE"
