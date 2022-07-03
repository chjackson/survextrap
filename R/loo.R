# Leave-one-out cross-validation for survextrap models
#
loo_survextrap <- function(x, standata){
    if (x$modelid != 1) stop("LOO not currently implemented for modelid 1")
    if (x$cure) stop("LOO not currently implemented for cure models")
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
    if (standata$ncovs>0){
        if (standata$nevent > 0) eta_event <- pars$loghr %*% t(standata$x_event)
        if (standata$nrcens > 0) eta_rcens <- pars$loghr %*% t(standata$x_rcens)
    } else eta_event <- eta_rcens <- 0
    if (standata$nevent > 0) {
        eta_event <- eta_event + standata$log_crude_event_rate +
            pars$gamma[,rep(1,standata$nevent)]
    }
    if (standata$nrcens > 0) {
        eta_rcens <- eta_rcens + standata$log_crude_event_rate +
            pars$gamma[,rep(1,standata$nrcens)]
    }
    eta_event <- unclass(eta_event)
    eta_rcens <- unclass(eta_rcens)
    ll_event <- log_dens(eta_event,  standata$basis_event,
                         pars$coefs, standata$cure, pars$cure_prob,
                         standata$ibasis_event, standata$modelid)
    ll_rcens <- log_surv(eta_rcens, standata$ibasis_rcens,
                         pars$coefs, standata$cure, pars$cure_prob,
                         standata$modelid)
    cbind(ll_event, ll_rcens)
}

log_dens <- function(eta, basis, coefs, cure, cure_prob, ibasis, modelid){
    log_haz(eta, basis, coefs, cure, cure_prob, ibasis, modelid) +
        log_surv(eta, ibasis, coefs, cure, cure_prob, modelid)
}

log_haz <- function(eta, basis, coefs, cure, cure_prob, ibasis, modelid){
    log(coefs %*% t(basis)) + eta
}

log_surv <- function(eta, ibasis, coefs, cure, cure_prob, modelid){
    - (coefs %*% t(ibasis)) * exp(eta)
}
