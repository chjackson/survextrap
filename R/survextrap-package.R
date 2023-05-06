#' The 'survextrap' package.
#'
#' @description For an introduction to and overview of the `survextrap` package, and full documentation, see
#'
#' \url{https://chjackson.github.io/survextrap/}
#' 
#' @docType package
#' @name survextrap-package
#' @useDynLib survextrap, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling nlist
#' @import ggplot2
#' @importFrom survival survfit
#' @importFrom posterior as_draws as_draws_matrix rhat ess_bulk
#' @importFrom gridExtra grid.arrange
#' @importFrom rstpm2 vuniroot
#' @importFrom splines2 mSpline iSpline
#' @importFrom survminer surv_summary
#' 
#' @importFrom stats as.formula delete.response dweibull integrate median model.matrix model.response optim predict pweibull quantile qbeta qgamma qnorm qt rlogis rnorm runif rbeta rgamma rt sd terms time var formula model.frame na.pass plogis qlogis get_all_vars reshape .getXlevels setNames
#'
NULL
