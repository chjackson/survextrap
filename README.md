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

* Extrapolations should not rely on conventional parametric forms (e.g. Weibull, log-normal, gamma...) that do not have interpretations as meaningful _mechanisms_ for how risk changes over time.

* Instead of selecting (or averaging) traditional parametric models, an _arbitrarily flexible_ parametric model should be used, that _adapts_ to give the optimal fit to the short-term and long-term data in combination.


### How it works 

* Bayesian multiparameter evidence synthesis is used to jointly model all sources of data and judgements 

* An M-spline is used to represent how the hazard changes through time (as in [rstanarm](https://arxiv.org/abs/2002.09633)).  The Bayesian fitting method automatically chooses the optimal level of smoothness and flexibility.  Spline "knots" should span the period covered by the data, and any period where there is a chance that the hazard may vary.

* A proportional hazards model is used to describe the relation of survival to predictors. 

* Mixture cure models are supported.

* It has an R interface, designed to be friendly to those familiar with standard R modelling functions.

* Stan is used under the surface to do MCMC (Hamiltonian Monte Carlo) sampling from the posterior distribution, in a similar fashion to [rstanarm](https://mc-stan.org/rstanarm/) and [survHE](https://CRAN.R-project.org/package=survHE). 

* Estimates and credible intervals for survival, hazard, mean and restricted mean survival can easily be extracted from the fitted model.


### Technical details of the methods

See `vignette("methods")`


### Examples of how to use it 

See `vignette("examples")`


## Development 

The package is in active development.  It can currently fit a large range of useful models, but it is not finished and is subject to be changed without warning.

Major things to do are:

* Relative survival / additive hazards models. 

* Thoughtful default priors for hazard changes without external data, e.g. in terms of orders of magnitude.

* Better thought-out knot choice, particularly with external data.

* More experience and examples of using it with real external data, including a vignette that lists how to implement all previously-suggested approaches for extrapolation with external data.

* Non-proportional hazards models.  This is expected to be computationally difficult, and I'm not sure of the best approach.

* Thorough testing, documentation and error handling.

If you want to try it out - feel free to install from github:

```{r}
remotes::install_github("chjackson/survextrap")
```

Please give feedback and suggestions if you do.  These can be posted on [github issues](https://github.com/chjackson/survextrap/issues), or [email](chris.jackson@mrc-bsu.cam.ac.uk).
