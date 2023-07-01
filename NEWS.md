# survextrap version 0.8.3 (2023/06/30)

* Fix for `mspline_spec` bug introduced in 0.8.2


# survextrap version 0.8.2 (2023/06/14)

* Individual data can now be excluded, and model fitted to external data alone


# survextrap version 0.8.1 (2023/05/25)

* Update to work with StanHeaders 2.26


# survextrap version 0.8 (2023/05/06)

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


# survextrap version 0.1

## Changes

See commit messages on GitHub for a history of developments to `survextrap`.
