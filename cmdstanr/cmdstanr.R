file <- file.path("inst", "stan", "survextrap.stan")
mod <- cmdstan_model(file)

## this makes a survextrap.exe file in inst/stan

## https://mc-stan.org/cmdstanr/articles/cmdstanr-internals.html#developing-using-cmdstanr
## To pre-compile all the models in a package, you may create top-level scripts
## configure and configure.win which run cmdstan_model() with compile = TRUE and
## save the compiled executables somewhere inside the inst/ folder of the
## package source. The instantiate package helps developers configure packages
## this way, and it documents other topics such as submitting to CRAN and
## administering CmdStan. Kevin Ushey’s configure package helps create and
## manage package configuration files in general.

## https://cran.r-project.org/web/packages/instantiate/index.html
## https://wlandau.github.io/instantiate/

## can we have both rstan and cmdstanr interfaces in the same package?
