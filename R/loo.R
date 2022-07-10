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
    pcure <- as.numeric(pars$pcure)
    ll_event <- log_dens(eta_event,  standata$basis_event,
                         pars$coefs, standata$cure, pcure,
                         standata$ibasis_event, standata$modelid)
    ll_rcens <- log_surv(eta_rcens, standata$ibasis_rcens,
                         pars$coefs, standata$cure, pcure,
                         standata$modelid)
    cbind(ll_event, ll_rcens)
}

log_dens <- function(eta, basis, coefs, cure, pcure, ibasis, modelid){
    log_haz(eta, basis, coefs, cure, pcure, ibasis, modelid) +
        log_surv(eta, ibasis, coefs, cure, pcure, modelid)
}

mspline_log_surv <- function(eta, ibasis, coefs) {
    - (coefs %*% t(ibasis)) * exp(eta);
}

mspline_log_haz <- function(eta, basis, coefs) {
    log(coefs %*% t(basis)) + eta
}

log_haz <- function(eta, basis, coefs, cure, pcure, ibasis, modelid){
    if (cure){
        if (modelid==1){
            base_logdens <- mspline_log_haz(eta, basis, coefs) + mspline_log_surv(eta, ibasis, coefs);
        } else  {
            base_logdens <- dweibull(basis, shape=coefs[,1], scale=exp(eta), log=TRUE)
        }
        res <- log(1 - pcure) + base_logdens -
            log_surv(eta, ibasis, coefs, cure, pcure, modelid)
    } else {
        if (modelid==1)
            res <- mspline_log_haz(eta, basis, coefs)
        else if (modelid==2)
            res <- dweibull(basis, shape=coefs[,1], scale=exp(eta), log=TRUE) -
                pweibull(ibasis, shape=coefs[,1], scale=exp(eta), log.p=TRUE, lower.tail=FALSE)
        else stop("unknown modelid")
    }
    res
}

log_surv <- function(eta, ibasis, coefs, cure, pcure, modelid){
    if (modelid==1)
        base_logsurv <- mspline_log_surv(eta, ibasis, coefs)
    else if (modelid==2)
        base_logsurv <- pweibull(ibasis, shape=coefs[,1], scale=exp(eta), log.p=TRUE, lower.tail=FALSE)
    else stop("unknown modelid")
    if (cure){
        res <- log(pcure + (1 - pcure)*exp(base_logsurv))
    } else {
        res <- base_logsurv
    }
    res
}
