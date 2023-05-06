
### SIMULATE DATA FOR CURE MODEL
set.seed(1)
n <- 200
x <- rep(c(0,1), each=round(n/2))
logor_cure <- 0.5
cure_prob <- plogis(qlogis(0.5) + logor_cure * x)
t <- numeric(n)
for (i in 1:n)
    t[i] <- flexsurvcure::rmixsurv(qweibull, n=1, theta=cure_prob[i], shape=1.5, scale=1.2)
censtime <- 10
curedata <- data.frame(t = pmin(t, censtime), status = as.numeric(t < censtime), x=x)
use_data(curedata, overwrite=TRUE)



### SIMULATE DATA FOR CURE MODEL: no covs on cure prob, covs on uncured
set.seed(1)
n <- 2000
x <- rep(c(0,1), each=round(n/2))
logor_cure <- 0
cure_prob <- plogis(qlogis(0.5) + logor_cure * x)
t <- numeric(n)
scale <- 1.2*exp(x)
for (i in 1:n)
  t[i] <- flexsurvcure::rmixsurv(flexsurv::qweibullPH, n=1, theta=cure_prob[i],
                                 shape=1.5, scale = scale[i])
censtime <- 10
curedata2 <- data.frame(t = pmin(t, censtime), status = as.numeric(t < censtime), x=x)
#use_data(curedata, overwrite=TRUE)
