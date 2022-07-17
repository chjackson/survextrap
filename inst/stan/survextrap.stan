/* Some of this code is from rstanarm (Sam Brilleman) */

functions {

    /**
    * Log hazard for M-spline model
    *
    * @param eta Vector, linear predictor
    * @param t Vector, event or censoring times
    * @param coefs Vector, M-spline coefficients
    * @return A vector
    */
    vector mspline_log_haz(vector eta, matrix basis, vector coefs) {
        return log(basis * coefs) + eta;
    }

    /**
    * Log survival and log CDF for M-spline model
    *
    * @param eta Vector, linear predictor
    * @param t Vector, event or censoring times
    * @param coefs Vector, M-spline coefficients
    * @return A vector
    */
    vector mspline_log_surv(vector eta, matrix ibasis, vector coefs) {
        vector[rows(eta)] res;
        res = - (ibasis * coefs) .* exp(eta);
        if (exp(res[1]) > 1) {
            reject("Probability > 1 computed. Not your fault - report a bug to the developer.");
        }
        return res;
    }

    vector mspline_log_dens(vector eta, matrix basis, matrix ibasis, vector coefs) {
        vector[rows(eta)] res;
        /* haz = dens / surv , loghaz = logdens - logsurv , logdens = loghaz + logsurv  */
        res = mspline_log_haz(eta, basis, coefs)  +
        mspline_log_surv(eta, ibasis, coefs);
        return res;
    }

    vector log_surv(vector eta, matrix ibasis, vector coefs,
		    data int cure, vector pcure, data int modelid) {
        vector[rows(eta)] res;
        vector[rows(eta)] base_logsurv;
        if (modelid==1){
            base_logsurv = mspline_log_surv(eta, ibasis, coefs);
        } else if (modelid==2) {
            for (i in 1:rows(eta)) {
                base_logsurv[i] = weibull_lccdf(ibasis[i,1] | coefs[1], exp(eta[i]));
            }
        }
        if (cure) {
            for (i in 1:rows(eta)) {
                res[i] = log(pcure[i] + (1 - pcure[i])*exp(base_logsurv[i]));
            }
        } else {
            res = base_logsurv;
        }
        return res;
    }

    vector log_haz(vector eta, matrix basis, vector coefs,
		   data int cure, vector pcure, matrix ibasis,
		   data int modelid) {
        vector[rows(eta)] res;
        vector[rows(eta)] base_logdens;
        vector[rows(eta)] base_loghaz;
        vector[rows(eta)] base_logsurv;
        if (cure) {
            if (modelid==1){
                base_logdens = mspline_log_dens(eta, basis, ibasis, coefs);
            } else {
                for (i in 1:rows(eta)){
                    base_logdens[i] = weibull_lpdf(basis[i,1] | coefs[1], exp(eta[i]));
                }
            }
            base_logsurv = log_surv(eta, ibasis, coefs, cure, pcure, modelid);
            for (i in 1:rows(eta)){
                res[i] = log(1 - pcure[i]) + base_logdens[i] - base_logsurv[i];
            }
        } else {
            if (modelid==1){
                base_loghaz = mspline_log_haz(eta, basis, coefs);
            } else {
                for (i in 1:rows(eta)){
                    base_loghaz[i] = weibull_lpdf(basis[i,1] | coefs[1], exp(eta[i])) -
                    weibull_lccdf(ibasis[i,1] | coefs[1], exp(eta[i]));
                }
            }
            res = base_loghaz;
        }
        return res;
    }

    vector log_dens(vector eta, matrix basis, vector coefs,
		    data int cure, vector pcure, matrix ibasis,
		    data int modelid){
        vector[rows(eta)] res;
        res = log_haz(eta, basis, coefs, cure, pcure, ibasis, modelid) +
        log_surv(eta, ibasis, coefs, cure, pcure, modelid);
        return res;
    }

    /**
    * Log-prior for intercept parameters
    *
    * @param gamma Real, the intercept parameter
    * @param dist Integer, the type of prior distribution
    * @param location Real, location of prior distribution
    * @param scale Real, scale for the prior distribution
    * @param df Real, df for the prior distribution
    * @return Nothing
    */
    real loghaz_lp(real gamma, int dist, real location, real scale, real df) {
        if (dist == 1)  // normal
        target += normal_lpdf(gamma | location, scale);
        else if (dist == 2)  // student_t
        target += student_t_lpdf(gamma | df, location, scale);
        /* else dist is 0 and nothing is added */
        return target();
    }

    real loghr_lp(vector loghr, int[] dist, vector location, vector scale, vector df) {
	for (i in 1:rows(loghr)){
	    if (dist[i] == 1)  // normal
		target += normal_lpdf(loghr[i] | location[i], scale[i]);
	    else if (dist[i] == 2)  // student_t
		target += student_t_lpdf(loghr[i] | df[i], location[i], scale[i]);
	}
        /* else dist is 0 and nothing is added */
        return target();
    }

}

data {
    int<lower=0> nevent;     // num. rows w/ an event (ie. not censored)
    int<lower=0> nrcens;     // num. rows w/ right censoring
    int<lower=0> nvars;      // num. aux parameters for baseline hazard
    int<lower=0> nextern;    // number of time points with external data
    int<lower=0> ncovs;      // number of covariate effects on the hazard
    int<lower=0> ncurecovs;      // number of covariate effects on the cure probability

    // log crude event rate / time (for centering linear predictor)
    real log_crude_event_rate;

    // basis matrices for M-splines / I-splines, without quadrature
    matrix[nevent,nvars] basis_event;  // at event time
    matrix[nevent,nvars] ibasis_event; // at event time
    matrix[nrcens,nvars] ibasis_rcens; // at right censoring time
    matrix[nextern,nvars] ibasis_ext_stop;  // at times with external data
    matrix[nextern,nvars] ibasis_ext_start; //
    matrix[nevent,ncovs] x_event; // matrix of covariate values
    matrix[nrcens,ncovs] x_rcens;
    matrix[nevent,ncurecovs] xcure_event; // matrix of covariate values on cure prob
    matrix[nrcens,ncurecovs] xcure_rcens;

    // external data describing knowledge about long-term survival
    // expressed as binomial outcomes of r survivors by t2 out of n people alive at t1
    int<lower=0> r_ext[nextern];
    int<lower=0> n_ext[nextern];
    matrix[nextern,ncovs] x_ext;
    matrix[nextern,ncurecovs] xcure_ext;

    vector[nvars-1] beta_mean; // logit of prior guess at basis weights (by default, those that give a constant hazard)
    int est_smooth;
    vector<lower=0>[1-est_smooth] smooth_sd_fixed;

    int cure;
    int prior_loghaz_dist;
    vector<lower=0>[3] prior_loghaz;
    vector<lower=0>[2] prior_cure;
    vector<lower=0>[2*est_smooth] prior_smooth;
    int<lower=0> prior_loghr_dist[ncovs];
    vector[ncovs] prior_loghr_location;
    vector[ncovs] prior_loghr_scale;
    vector[ncovs] prior_loghr_df;
    int<lower=0> prior_logor_cure_dist[ncurecovs];
    vector[ncurecovs] prior_logor_cure_location;
    vector[ncurecovs] prior_logor_cure_scale;
    vector[ncurecovs] prior_logor_cure_df;

    int modelid;
}

parameters {
    real gamma[1];
    vector[ncovs] loghr;
    vector[nvars-1] beta_err;
    vector<lower=0>[est_smooth] smooth_sd;
    vector<lower=0,upper=1>[cure] pcure;
    vector[ncurecovs] logor_cure;
}


transformed parameters {
    vector[nvars] beta;
    vector[nvars] coefs; // constrained coefs for M-splines

    if (est_smooth)
    beta = append_row(0, beta_mean + beta_err*smooth_sd[1]);
    else
    beta = append_row(0, beta_mean + beta_err*smooth_sd_fixed[1]);
    if (modelid==1){
        coefs = softmax(beta);
    } else {
        coefs[1] = exp(beta_err[1]);
        coefs[2] = 0; // dummy, unused
    }

}

model {
    vector[nevent] eta_event; // for events
    vector[nrcens] eta_rcens; // for right censored
    vector[nextern] eta_extern; // for external data
    real dummy;
    real cp;
    vector[nextern] p_ext_stop; // unconditional survival prob at external time points
    vector[nextern] p_ext_start; //
    vector[nevent] pcure_event; //
    vector[nrcens] pcure_rcens; //
    vector[nextern] pcure_extern; //

    if (ncovs > 0) {
        // does x * beta,    matrix[n,K] * vector[K] is this a matrix product
        if (nevent > 0) eta_event = x_event * loghr;
        if (nrcens > 0) eta_rcens = x_rcens * loghr;
        if (nextern > 0) eta_extern = x_ext * loghr;
    } else {
        if (nevent > 0) eta_event = rep_vector(0.0, nevent);
        if (nrcens > 0) eta_rcens = rep_vector(0.0, nrcens);
        if (nextern > 0) eta_extern = rep_vector(0.0, nextern);
    }

    // add on log crude event rate / time (helps to center intercept)
    if (nevent > 0) eta_event += log_crude_event_rate;
    if (nrcens > 0) eta_rcens += log_crude_event_rate;
    if (nextern > 0) eta_extern += log_crude_event_rate;

    // add on intercept to linear predictor
    if (nevent > 0) eta_event += gamma[1];
    if (nrcens > 0) eta_rcens += gamma[1];
    if (nextern > 0) eta_extern += gamma[1];

    if (cure) cp = pcure[1]; else cp = 0;
    pcure_event = rep_vector(cp, nevent);
    pcure_rcens = rep_vector(cp, nrcens);
    pcure_extern = rep_vector(cp, nextern);

    if (ncurecovs > 0){
        if (nevent > 0) pcure_event = inv_logit(logit(pcure_event) + xcure_event * logor_cure);
        if (nrcens > 0) pcure_rcens = inv_logit(logit(pcure_rcens) + xcure_rcens * logor_cure);
        if (nextern > 0) pcure_extern = inv_logit(logit(pcure_extern) + xcure_ext * logor_cure);
    }

    if (nevent > 0) target +=  log_dens(eta_event,  basis_event, coefs, cure, pcure_event, ibasis_event, modelid);
    if (nrcens > 0) target +=  log_surv(eta_rcens, ibasis_rcens, coefs, cure, pcure_rcens, modelid);

    if (nextern > 0) {
        p_ext_stop = exp(log_surv(eta_extern, ibasis_ext_stop, coefs, cure, pcure_extern, modelid));
        p_ext_start = exp(log_surv(eta_extern, ibasis_ext_start, coefs, cure, pcure_extern, modelid));
        target += binomial_lpmf(r_ext | n_ext, p_ext_stop ./ p_ext_start);
    }

    // log prior for baseline log hazard
    dummy = loghaz_lp(gamma[1], prior_loghaz_dist, prior_loghaz[1],
		      prior_loghaz[2], prior_loghaz[3]);

    // log prior for covariates on the log hazard
    dummy = loghr_lp(loghr, prior_loghr_dist, prior_loghr_location,
        		     prior_loghr_scale, prior_loghr_df);

    // prior for spline coefficient random effect term
    if (modelid==1){
        beta_err ~ logistic(0, 1);
    }

    // prior for cure fraction
    if (cure) {
        pcure ~ beta(prior_cure[1], prior_cure[2]);
    }
    if (ncurecovs > 0){
	// log prior for covariates on the log odds of cure
	dummy = loghr_lp(logor_cure, prior_logor_cure_dist, prior_logor_cure_location,
			 prior_logor_cure_scale, prior_logor_cure_df);
    }
    
    if (est_smooth){
        smooth_sd ~ gamma(prior_smooth[1], prior_smooth[2]);
    }
}

generated quantities {
    // transformed intercept.  TODO Why not call this eta?  Or just loghaz cos it is log of eta in the methods vig
    real alpha = log_crude_event_rate + gamma[1];
}
