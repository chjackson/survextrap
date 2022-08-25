#' The 'survextrap' package.
#'
#' @description For an introduction to and an overview of the `survextrap` package, see
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
#' @importFrom stats as.formula delete.response dweibull integrate median model.matrix model.response optim predict pweibull quantile reformulate rlogis rnorm runif sd terms time var formula model.frame na.pass plogis qlogis get_all_vars reshape
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.26.11. https://mc-stan.org
#'
NULL
