# Version 0.8.15 (2024/06/29)

* Random walk prior for the coefficients changed to a first-order weighted model, following [Phillippo et al.](https://arxiv.org/abs/2401.12640), so that the expected smoothness accounts for the knot spacing.


# Version 0.8.14 (2024/06/13)

* Added discount rate argument `disc_rate` to `mean` and `rmst()` functions (PR #12)


# Version 0.8.13 (2024/03/05)

* Bug fix for models with external data and either nonproportional hazards or covariates on cure probability. 


# Version 0.8.12 (2024/01/06)

* `irmst`, `hazard_ratio` and `hrtime` tidied.   Fix of bug in `irmst` introduced in v0.8.9, and all now have a `sample` argument.


# Version 0.8.11 (2023/12/28)

* Outputs standardised over the covariate values of a reference population can now be produced using the new `standardise_to()` function.


# Version 0.8.10 (2023/12/24)

* `loo` no longer ignores covariates on the cure probability in cure models.

* More internal refactoring of code to produce outputs.


# Version 0.8.9 (2023/12/22)

* Cleaned and simplified internal code for model outputs (`rmst`, `survival` etc.). 

* `rmst` now outputs a tibble, and credible limits in the default output are renamed to `lower` and `upper` for consistency with `survival` and `hazard`.

* Arrays now allowed for `coef` argument in `dsurvspline` and related functions.


# Version 0.8.8 (2023/12/19)

* Background hazards in additive hazards models can now be stratified, using the new `backhaz_strata` argument to `survextrap()`.  This argument names stratifying variables (e.g. age group) that should be in both the background hazards and the individual and/or external data.

* Bug fix for `hsurvmspline_wane` and `dsurvmspline_wane` where offsets were not implemented correctly.  This affected hazard predictions for models that included both background hazards and treatment effect waning.


# Version 0.8.7 (2023/10/29)

* A random walk prior can now be specified for the spline coefficients, by `smooth_model = "random_walk"`. 


# Version 0.8.6 (2023/10/24)

* Non-proportional hazards models can now be applied to a subset of the covariates, by supplying a formula as the `nonprop` argument to `survextrap`. 

* Prior simulation functions such as `prior_sample` now require an explicit covariate model, specified through `formula` and `newdata`, rather than a design matrix `X`.  These functions now fully support nonproportional hazards models and cure probabilities with associated regression models.


# Version 0.8.5 (2023/09/30)

* Code for constructing `mspline` specification objects tidied, with new `mspline_init` and `mspline_list_init` functions.


# Version 0.8.4 (2023/09/09)

* Deprecated Stan array syntax updated.  The package now requires `rstan` version 2.26.


# Version 0.8.3 (2023/06/30)

* Fix for `mspline_spec` bug introduced in 0.8.2.


# Version 0.8.2 (2023/06/14)

* Individual data can now be excluded, so that the model is fitted to external data alone.


# Version 0.8.1 (2023/05/25)

* Update to work with StanHeaders 2.26.


# Version 0.8 (2023/05/06)

The [GitHub commit](https://github.com/chjackson/survextrap/commit/1668f40604d9dc62d83a698c735275506474aa03) on 6 May 2023 implemented several changes that are listed here.  `survextrap` is now in "beta" status.  All major features are now implemented, but there may still be some lack of polish.


## New features 

* A new vignette `cetuximab`, giving an in-depth case study of a realistic use of `survextrap`.  To be included in a forthcoming paper.

* Improvements to usability and stability of M-splines:

    - `iknots` and `bknots` are replaced by a single `knots`, and the lower boundary knot is always fixed at zero.
	
	- `add_knots` can be used to add user-defined knots to the default ones.  See the case study.

    - New function `mspline_spec` to define M-spline knots based on data, in advance of fitting any models. 

    - M-splines can now be smoother at the final knot, through an option `bsmooth`, which is on by default.  Thanks to Iain Timmins for the suggestion.

* Incremental restricted mean survival time (`irmst`) added.

* Any summary of the posterior distribution can be produced by supplying a list of summary functions (e.g. `mean`, `quantiles`) to output functions (such as `rmst` and `survival`).

* Leave-one-out cross validation for external data.


## Bug fixes and improvements

* Bug fix for `hsurvmspline()` (hence `hazard()`) which was using the wrong tail end of `psurvmspline()`.

* Waning models no longer ignore non-proportional hazards and cure.  Waning vignette merged into the methods vignette.

* Leave-one-out cross validation now supports models with background hazards.

* More consistent naming of arguments and variables (e.g. `hscale` for hazard scale, and `hsd` for hazard smoothing parameter).

* More information given when printing fitted model objects.

* Thorough documentation polishing and code tidying.


# Version 0.1

## Changes

See commit messages on GitHub for a history of developments to `survextrap`.
