# survextrap 

`survextrap` is an R package under development, to model survival from a combination of 

1. A standard individual-level, right-censored survival dataset

2. "External" data sources in the following aggregate "count" form:

<table> 
<tr>
<th colspan="2">Follow-up period </th>
<th colspan="2">Number</th>
</tr> 
<tr><th>Start time $t$</th><th>End time $u$</th><th>Alive at $t$</th><th>Still alive at $u$</th></tr>

<tr>
<td> $t_{1}$ </td>
<td> $u_{1}$ </td>
<td> $n_{1}$ </td>
<td> $r_{1}$ </td>
</tr>

<tr>
<td> $t_{2}$ </td>
<td> $u_{2}$ </td>
<td> $n_{2}$ </td>
<td> $r_{2}$ </td>
</tr>

<tr>
<td>etc...</td>
<td></td>
<td></td>
<td></td>

</tr>

</table>

Any number of rows can be supplied for the "external" data, and the time intervals do not have to be distinct or exhaustive. 

The package has been developed under the expectation that many forms of external data that might be useful for survival extrapolation (such as population data, registry data or elicited judgements) can be manipulated into this common "count" form.

### Principles

* Extrapolations from short-term individual level data should be done using _explicit data or judgements_ about how risk will change over time. 

* Extrapolations should not rely on standard parametric forms (e.g. Weibull, log-normal, gamma...) that are only used out of convention and do not have interpretations as plausible _mechanisms_ for how risk will change over time.

* Instead of selecting (or averaging) traditional parametric models, an _arbitrarily flexible_ parametric model should be used, that _adapts_ to give the optimal fit to the short-term and long-term data in combination.


### How it works 

* Bayesian multiparameter evidence synthesis is used to jointly model all sources of data and judgements.

* An M-spline is used to represent how the hazard changes through time (as in [rstanarm](https://arxiv.org/abs/2002.09633)).  The Bayesian fitting method automatically chooses the optimal level of smoothness and flexibility.  Spline "knots" should span the period covered by the data, and any future period where there is a chance that the hazard may vary.  Then if there is no data in the future period, the uncertainty will be acknowledged and the predicted hazards will have wide credible intervals.

* A proportional hazards model or a flexible non-proportional hazards model can be used to describe the relation of survival to predictors. 

* Mixture cure and relative survival models are supported.

* It has an R interface, designed to be friendly to those familiar with standard R modelling functions.

* Stan is used under the surface to do MCMC (Hamiltonian Monte Carlo) sampling from the posterior distribution, in a similar fashion to [rstanarm](https://mc-stan.org/rstanarm/) and [survHE](https://CRAN.R-project.org/package=survHE). 

* Estimates and posterior credible intervals / samples for survival, hazard, mean and restricted mean survival can easily be extracted.


### Technical details of the methods

See `vignette("methods")`


### Examples of how to use it 

See `vignette("examples")`


## Development 

The package is in active development.  It can currently fit a large range of useful models, but it is not finished and is subject to be changed without warning.

Major things to do are:

* Thoughtful default priors for hazard changes without external data, e.g. in terms of orders of magnitude.  Better thought-out knot choice, particularly with external data.  Empirical work to show the impact of priors and knot choice.

* More experience and examples of using it with real external data, including a vignette that lists how to implement all previously-suggested approaches for extrapolation with external data.

* Recommendations for how to use cure and relative survival together, and consider whether these need to be extended.  For example, uncertainty in the background hazard, with non-proportional hazards with study population (additive or more flexible?)

* Explain relation to dynamic GLMs as used by [Kearns et al](https://doi.org/10.1177%2F0272989X19873661) and related papers. Can "random walk" principles used in their dynamic model inspire a sensible prior for hazard changes?

* Thorough testing, documentation and error handling.

If you want to try it out - feel free to install the development version as: 

```{r}
install.packages("survextrap", repos=c('https://chjackson.r-universe.dev', 'https://cloud.r-project.org'))
```

Please give feedback and suggestions if you do.  These can be posted on [github issues](https://github.com/chjackson/survextrap/issues), or [email](mailto:chris.jackson@mrc-bsu.cam.ac.uk).

<!-- badges: start -->
[![R-CMD-check](https://github.com/chjackson/survextrap/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/chjackson/survextrap/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/chjackson/survextrap/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/chjackson/survextrap/actions/workflows/test-coverage.yaml)
<!-- badges: end -->
