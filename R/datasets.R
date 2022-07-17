##' Colon cancer survival data
##' 
##' Survival times of 185 patients with primary node positive breast
##' cancer.  This is the data provided in the `survival` package as
##' `colon`, but this is restricted to outcome of recurrence-free
##' survival, artificially censored at 3 years, and taking a 20%
##' random sample.
##' 
##' @format Documented in \code{\link[survival]{colon}}. 
##'
##' @references See \code{\link[survival]{colon}}. 
##'
##' @source See \code{\link[survival]{colon}}. 
##'
##' @keywords datasets
"colons"



##' Simulated data for testing cure models
##' 
##' Survival times of 200 fake people.
##' The cure probability is 0.5 for `x=0`, and the log
##' odds ratio for cure is 0.5, so that the cure probability is about 0.62 for `x=1`.
##' The uncured population have survival times distributed as a Weibull with shape 1.5 and scale 1.2.
##' 
##' @format
##'
##' `t` Survival times, right-censored at 10 years. 
##' 
##' `status`. 1 for an observed death and 0 for censoring.
##'
##' `x`. A numeric variable with value of 0 for 100 individuals, and 1 for the other 100. 
##'
##' @source Simulated. 
##'
##' @keywords datasets
"curedata"
