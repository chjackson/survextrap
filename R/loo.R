# Leave-one-out cross-validation for survextrap models
#
loo_survextrap <- function(x, standata){
    ll <- loglik_ipd(x, standata)
    chains <- x$stanfit@sim$chains
    iter <- NROW(ll)
    chain_id <- rep(1:chains, each=iter/chains)
    r_eff <- loo::relative_eff(exp(ll), chain_id=chain_id, cores=1)
    res <- loo::loo(ll, r_eff = r_eff)
    res
}

# Mimic the Stan code

loglik_ipd <- function(x, standata){
    pars <- get_pars(x, newdata=NULL)
    coefs_event <- get_coefs_bycovs(x, get_draws(x), X=standata$x_event)
    coefs_rcens <- get_coefs_bycovs(x, get_draws(x), X=standata$x_rcens)

    if (standata$ncovs>0){
        if (standata$nevent > 0) alpha_event <- pars$loghr %*% t(standata$x_event)
        if (standata$nrcens > 0) alpha_rcens <- pars$loghr %*% t(standata$x_rcens)
    } else alpha_event <- alpha_rcens <- 0
    if (standata$nevent > 0) {
        alpha_event <- alpha_event + standata$prior_loghaz[1] +
            pars$gamma[,rep(1,standata$nevent)]
    }
    if (standata$nrcens > 0) {
        alpha_rcens <- alpha_rcens + standata$prior_loghaz[1] +
            pars$gamma[,rep(1,standata$nrcens)]
    }
    alpha_event <- unclass(alpha_event)
    alpha_rcens <- unclass(alpha_rcens)
    pcure <- as.numeric(pars$pcure)
    ll_event <- log_dens(alpha_event,  standata$basis_event,
                         coefs_event, standata$cure, pcure,
                         standata$ibasis_event, modelid=1)
    ll_rcens <- log_surv(alpha_rcens, standata$ibasis_rcens,
                         coefs_rcens, standata$cure, pcure,
                         modelid=1)
    cbind(ll_event, ll_rcens)
}

log_dens <- function(alpha, basis, coefs, cure, pcure, ibasis, modelid){
    log_haz(alpha, basis, coefs, cure, pcure, ibasis, modelid) +
        log_surv(alpha, ibasis, coefs, cure, pcure, modelid)
}

mspline_log_surv <- function(alpha,
                             ibasis, # nobs x df
                             coefs # niter x nobs x df
                             ) {
    ibasis_arr <- array(ibasis, dim=c(1, dim(ibasis)))[rep(1,dim(coefs)[1]),,] # niter x nobs x df
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

log_haz <- function(alpha, basis, coefs, cure, pcure, ibasis, modelid){
  if (cure){
    base_logdens <- mspline_log_haz(alpha, basis, coefs) + mspline_log_surv(alpha, ibasis, coefs);
    res <- log(1 - pcure) + base_logdens -
      log_surv(alpha, ibasis, coefs, cure, pcure, modelid)
  } else {
    res <- mspline_log_haz(alpha, basis, coefs)
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
