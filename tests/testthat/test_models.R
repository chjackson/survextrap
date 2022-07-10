### TODO
### Use a dataset in the survival package, e.g. why not colon or ovarian
### Not full MCMC.  Optimisation will do
### Though need at least one test of postest with MCMC and VB
### Use small samples to run mean, RMST, survival and hazard after each one.
### postest_test() function for all if they will look the same

### Purpose of this is to test postest with different model forms
### Neednt do it for everything, only where necessary
### Could be additional tests for features of output functions
postest_test <- function(x, newdata=NULL){
    niter <- 10
    t <- c(1, 5)
    summary(x)
    mean(x, newdata=newdata, niter=niter)
    rmst(x, t=t, newdata=newdata, niter=niter)
    survival(x, newdata=newdata, times=t)
    hazard(x, newdata=newdata, times=t)
    invisible()
}
library(dplyr)
colonc <- colon |>
    mutate(years = time/365.25) |>
    filter(etype == 1) |>
    mutate(status = ifelse(years > 3, 0, status),
           years = ifelse(years > 3, 3, years))


### Basic spline model, no covariates
mod <- survextrap(Surv(years, status) ~ 1, data=colonc, fit_method="opt")
postest_test(mod)

### Basic spline model, with covariates
mod <- survextrap(Surv(years, status) ~ rx, data=colonc, fit_method="opt")
nd <- data.frame(rx = c("Obs", "Lev+5FU"))
postest_test(mod, nd)

### Weibull model, no covariates
mod <- survextrap(Surv(years, status) ~ 1, data=colonc, fit_method="opt", modelid="weibull")
postest_test(mod)

### Weibull model, covariates
mod <- survextrap(Surv(years, status) ~ rx, data=colonc, fit_method="opt", modelid="weibull")
postest_test(mod, nd)
# flexsurv::flexsurvreg(Surv(years, status) ~ rx, data=colonc, dist="weibull")

## Weibull, MCMC fit
mod <- survextrap(Surv(years, status) ~ 1, data=colonc, fit_method="mcmc",
                  modelid="weibull", chains=1, iter=1000)
postest_test(mod, nd)


### Changing various settings in basic model
### Internal knots, boundary knots
### Spline prior mean
### EB for smooth variance
### Fixed smooth variance

### With external data

### Fitting method
### VB

### SIMULATE DATA FOR CURE MODEL
### TODO remove dependence, since that sucks in flexsurvcure, or do in Suggests.
### It's the imports we want to keep slim
set.seed(1)
x <- rep(c(0,1), each=100)
logor_cure <- 0.5
cure_prob <- plogis(qlogis(0.5) + logor_cure * x)
t <- numeric(200)
for (i in 1:200)
    t[i] <- flexsurvcure::rmixsurv(qweibull, n=1, theta=cure_prob[i], shape=1.5, scale=1.2)
censtime <- 10
cure_df <- data.frame(t = pmin(t, censtime), status = as.numeric(t < censtime), x=x)
nd <- data.frame(x=c(0,1))

### Cure model, no covariates on anything
cmod0 <- survextrap(Surv(t, status) ~ 1, data=cure_df, chains=1, cure=TRUE, iter=1000, fit_method="opt")
postest_test(cmod0)

### Cure model, covariates on noncured model
cmod1 <- survextrap(Surv(t, status) ~ x, data=cure_df, chains=1, cure=TRUE, iter=1000, fit_method="opt")
postest_test(cmod1, newdata=nd)

### Cure model, covariates on cured fraction
cmod2 <- survextrap(Surv(t, status) ~ 1, data=cure_df, chains=1, cure=~x, iter=1000, fit_method="opt")
postest_test(cmod2, newdata=nd)
