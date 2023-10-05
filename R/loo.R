## Leave-one-out cross-validation for survextrap models
## Implemented by mimicking the Stan code for likelihood computation in R 

loo_survextrap <- function(x, standata, loglik_fn){
    ll <- loglik_fn(x, standata)
    chains <- x$stanfit@sim$chains
    iter <- NROW(ll)
    chain_id <- rep(1:chains, each=iter/chains)
    r_eff <- loo::relative_eff(exp(ll), chain_id=chain_id, cores=1)
    res <- loo::loo(ll, r_eff = r_eff)
    res
}

## Log-likelihood for the individual-level data

loglik_ipd <- function(x, standata){
    pars <- get_pars(x, newdata=NULL)
    coefs_event <- get_coefs_bycovs(x, get_draws(x), X=standata$x_event)
    coefs_rcens <- get_coefs_bycovs(x, get_draws(x), X=standata$x_rcens)

    if (standata$ncovs>0){
        if (standata$nevent > 0) alpha_event <- pars$loghr %*% t(standata$x_event)
        if (standata$nrcens > 0) alpha_rcens <- pars$loghr %*% t(standata$x_rcens)
    } else alpha_event <- alpha_rcens <- 0
    if (standata$nevent > 0) {
        alpha_event <- alpha_event + standata$prior_hscale[1] +
            pars$gamma[,rep(1,standata$nevent)]
    }
    if (standata$nrcens > 0) {
        alpha_rcens <- alpha_rcens + standata$prior_hscale[1] +
            pars$gamma[,rep(1,standata$nrcens)]
    }
    alpha_event <- unclass(alpha_event)
    alpha_rcens <- unclass(alpha_rcens)
    pcure <- as.numeric(pars$pcure)
    ll_event <- log_dens(alpha_event,  standata$basis_event,
                         coefs_event, standata$cure, pcure,
                         standata$ibasis_event, modelid=1,
                         standata$relative,
                         standata$backhaz_event)
    ll_rcens <- log_surv(alpha_rcens, standata$ibasis_rcens,
                         coefs_rcens, standata$cure, pcure,
                         modelid=1)
    cbind(ll_event, ll_rcens)
}

## Aggregated likelihood is r ~ Bin(n, pstop/pstart)
## For LOO, express at the individual level as
## ri ~ Bern((pstop/pstart)^ri (1 - pstop/pstart)^(1-ri))
## for each 0/1 event in the disaggregated data

loglik_external <- function(x, standata){
  pars <- get_pars(x, newdata=NULL)
  coefs <- get_coefs_bycovs(x, get_draws(x), X=standata$x_ext)
  if (standata$ncovs>0){
    alpha <- pars$loghr %*% t(standata$x_ext)
  } else alpha <-0
  alpha <- alpha + standata$prior_hscale[1] + pars$gamma[,rep(1,standata$nextern)]
  alpha <- unclass(alpha)
  pcure <- as.numeric(pars$pcure)
  y <- rep(rep(c(0, 1), each=standata$nextern),
           c(standata$n_ext - standata$r_ext, standata$r_ext))
  inds <- rep(1:standata$nextern, standata$n_ext)
  niter <- attr(pars, "niter")
  lp_stop <- log_surv(alpha,  standata$ibasis_ext_stop, coefs,
                      standata$cure, pcure, modelid=1) +
    rep(log(standata$backsurv_ext_stop), each=niter)
  lp_start <- log_surv(alpha,  standata$ibasis_ext_start, coefs,
                       standata$cure, pcure, modelid=1) +
    rep(log(standata$backsurv_ext_start), each=niter)
  lp_stop <- lp_stop[,inds] # niter rows x ni cols
  lp_start <- lp_start[,inds]
  y <- matrix(y, nrow=1)[rep(1,niter),] # niter x ni
  y*(lp_stop - lp_start) + (1 - y)*log(1 - exp(lp_stop - lp_start))
}

log_dens <- function(alpha, basis, coefs, cure, pcure, ibasis, modelid,
                     relative, backhaz){
  log_haz(alpha, basis, coefs, cure, pcure, ibasis, modelid,
          relative, backhaz) +
    log_surv(alpha, ibasis, coefs, cure, pcure, modelid)
}

mspline_log_surv <- function(alpha,
                             ibasis, # nobs x df
                             coefs # niter x nobs x df
                             ) {
    ibasis_arr <- array(ibasis, dim=c(1, dim(ibasis)))[rep(1,dim(coefs)[1]),,,drop=FALSE] # niter x nobs x df
    coefs_x_ibasis <- apply(coefs * ibasis_arr, c(1,2), sum)
    - (coefs_x_ibasis) * exp(alpha) # niter x nobs
}

mspline_log_haz <- function(alpha,
                            basis, # nobs x df
                            coefs # niter x nobs x df
                            ) {
    basis_arr <- array(basis, dim=c(1, dim(basis)))[rep(1,dim(coefs)[1]),,] # niter x nobs x df
    coefs_x_basis <- apply(coefs * basis_arr, c(1,2), sum)
    log(coefs_x_basis) + alpha # niter x nobs
}

log_haz <- function(alpha, basis, coefs, cure, pcure, ibasis, modelid,
                    relative, backhaz){
  if (cure){
    base_logdens <- mspline_log_haz(alpha, basis, coefs) + mspline_log_surv(alpha, ibasis, coefs);
    res <- log(1 - pcure) + base_logdens -
      log_surv(alpha, ibasis, coefs, cure, pcure, modelid)
  } else {
    res <- mspline_log_haz(alpha, basis, coefs)
  }
  if (relative){
    niter <- nrow(res)
    nobs <- ncol(res)
    backhaz <- matrix(backhaz, nrow=niter, ncol=nobs, byrow=TRUE)
    res <- log(backhaz + exp(res))
  }
  res
}

log_surv <- function(alpha, ibasis, coefs, cure, pcure, modelid){
  base_logsurv <- mspline_log_surv(alpha, ibasis, coefs)
  if (cure){
    res <- log(pcure + (1 - pcure)*exp(base_logsurv))
  } else {
    res <- base_logsurv
  }
  res
}
