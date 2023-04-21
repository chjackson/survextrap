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


##' Datasets for evaluation of cetuximab in head and neck cancer
##'
##' Datasets for evaluation of cetuximab in head and neck cancer, as previously
##' analysed by Guyot et al. (2017) to demonstrate models for survival extrapolation
##' with Bayesian evidence synthesis.
##'
##' @aliases cetux cetux_seer cetux_bh
##'
##' 
##' @format \code{cetux} contains synthetic individual-level survival
##'   data, generated (using the method of Guyot et al 2012) to be
##'   consistent with the Kaplan-Meier estimates of survival published
##'   by Bonner et al. (2006).
##' Columns \code{months}, \code{years}, give survival in months or
##' years since the date of diagnosis of head and neck cancer.
##' \code{d} is a numeric vector where 0 indicates censoring, and 1
##' indicates death at this time.  \code{treat} is a factor indicating
##' the treatment group (Cetuximab or Control); both groups also
##' received radiotherapy.
##'
##' \code{cetux_seer} Estimates of conditional survival from registry data, 
##' matched to the Bonner trial population by age, gender, cancer site, and date of diagnosis.
##' From the "Surveillance Epidemiology and End Results" (SEER) database.
##' Each line gives counts of \code{r} survivors up to \code{stop} years, 
##' given \code{n} people alive at \code{start}.  \code{haz} is the
##' corresponding constand hazard estimate over this period. 
##'
##' \code{cetux_bh} Mortality rates for the population of the
##' USA, matched by age and sex to the patients from the Bonner trial.
##' 80% are male, and the median age is 57 (range 34 to 83).  Hence
##' the \eqn{i}th row is a weighted average of the male and 
##' and female mortality rates for age \eqn{57 + i - 1}.
##'
##' See Guyot et al. (2017) for more details of each of these.
##'
##' @references 
##'
##' Guyot, P., Ades, A.E., Beasley, M., Lueza, B., Pignon, J.P. and
##' Welton, N.J., (2017). Extrapolation of survival curves from cancer
##' trials using external information. Medical Decision Making, 37(4),
##' pp.353-366.
##'
##' Bonner, J.A., Harari, P.M., Giralt, J., Azarnia, N., Shin, D.M.,
##' Cohen, R.B., Jones, C.U., Sur, R., Raben, D., Jassem, J. and Ove,
##' R., (2006) Radiotherapy plus cetuximab for squamous-cell carcinoma
##' of the head and neck. New England Journal of Medicine, 354(6),
##' pp.567-578.
##'
##' Guyot, P., Ades, A.E., Ouwens, M.J. and Welton, N.J.,
##' (2012) Enhanced secondary analysis of survival data: reconstructing
##' the data from published Kaplan-Meier survival curves. BMC Medical
##' Research Methodology, 12, pp.1-13.
##'
##' @name cetux
NULL

##' @rdname cetux
"cetux"

##' @rdname cetux
"cetux_seer"

##' @rdname cetux
"cetux_bh"
