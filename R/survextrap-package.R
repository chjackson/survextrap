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
#' @import dplyr
#' @import ggplot2
#' @importFrom survival survfit
#' @importFrom tidyr extract pivot_longer
#' @importFrom posterior as_draws
#' @importFrom magrittr "%>%"
#' @importFrom gridExtra grid.arrange
#' @importFrom rstpm2 vuniroot
#' @importFrom splines2 mSpline iSpline
#' @importFrom survminer surv_summary
#' @importFrom tidybayes gather_rvars median_qi
#' @importFrom stats as.formula delete.response dweibull integrate median model.matrix model.response optim predict pweibull quantile reformulate rlogis rnorm runif sd terms time var formula model.frame na.pass plogis qlogis get_all_vars
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.26.11. https://mc-stan.org
#'
NULL
