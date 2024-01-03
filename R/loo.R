#' Leave-one-out cross-validation for survextrap models Implemented by
#' mimicking the Stan code for likelihood computation in R
#'
#' @param x Result of call to \code{survextrap}
#'
#' @param standata List of data supplied in the call to \code{rstan}'s
#'   MCMC or other fitting function
#'
#' @param loglik_fn Either \code{loglik_ipd} or \code{loglik_external}
#'
#' @return Result of call to \code{loo::loo}
#'
#' @noRd 
loo_survextrap <- function(x, standata, loglik_fn){
  ll <- loglik_fn(x, standata)
  chains <- x$stanfit@sim$chains
  iter <- NROW(ll)
  chain_id <- rep(1:chains, each=iter/chains)
  r_eff <- loo::relative_eff(exp(ll), chain_id=chain_id, cores=1)
  res <- loo::loo(ll, r_eff = r_eff)
  res
}

#' Posterior log-likelihoods for the individual-level data
#'
#' @inheritParams loo_survextrap
#'
#' @return Matrix with one column per individual-level observation
#'   (ordered with event times before censoring times) and one row per
#'   MCMC sample.
#'
#' @noRd 
loglik_ipd <- function(x, standata){
  stanmat <- get_draws(x)
  if (standata$nevent > 0) {
    alpha_event <- get_alpha_bycovs(x, stanmat, X=standata$x_event)
    coefs_event <- get_coefs_bycovs(x, stanmat, X=standata$xnph_event)
    pcure_event <- get_pcure_bycovs(x, stanmat, X=standata$x_event)
    ll_event <- log_dens(alpha_event,  standata$basis_event,
                         coefs_event, standata$cure, pcure_event,
                         standata$ibasis_event, modelid=1,
                         standata$relative,
                         standata$backhaz_event)
  } else ll_event <- NULL
  if (standata$nrcens > 0) {
    alpha_rcens <- get_alpha_bycovs(x, stanmat, X=standata$x_rcens)
    coefs_rcens <- get_coefs_bycovs(x, stanmat, X=standata$xnph_rcens)
    pcure_rcens <- get_pcure_bycovs(x, stanmat, X=standata$x_rcens)
    ll_rcens <- log_surv(alpha_rcens, standata$ibasis_rcens,
                         coefs_rcens, standata$cure, pcure_rcens,
                         modelid=1)
  } else ll_rcens <- NULL
  cbind(ll_event, ll_rcens)
}

#' Posterior log-likelihoods for the external data
#' 
#' @details
#'
#' Aggregated likelihood is \eqn{r \sim Bin(n, pstop/pstart)}.
#'
#' For LOO, this is expressed at the individual level, as
#'
#' \deqn{ri \sim Bern((pstop/pstart)^ri (1 - pstop/pstart)^(1-ri))}
#'
#' for each 0/1 event in the disaggregated data
#'
#' @inheritParams loo_survextrap
#'
#' @return Matrix with one column per individual in the disaggregated
#'   external data, and one row per MCMC sample.
#'
#' @noRd 
loglik_external <- function(x, standata){
  stanmat <- get_draws(x)
  coefs <- get_coefs_bycovs(x, stanmat, X=standata$xnph_ext)
  alpha <- get_alpha_bycovs(x, stanmat, X=standata$x_ext)
  pcure <- get_pcure_bycovs(x, stanmat, X=standata$x_ext)
  y <- rep(rep(c(0, 1), each=standata$nextern),
           c(standata$n_ext - standata$r_ext, standata$r_ext))
  inds <- rep(1:standata$nextern, standata$n_ext)
  niter <- nrow(stanmat)
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
