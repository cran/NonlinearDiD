# NonlinearDiD 0.2.0

* Added support for repeated cross-section staggered DiD designs.
* Added `data_type = "panel"` and `data_type = "repeated_cross_section"`
  options.
* Made `idname` optional for repeated cross-section designs.
* Added support for sampling weights through `weightsname`.
* Added support for clustered inference through `cluster_var`.
* Added examples for binary repeated cross-section outcomes with staggered
  treatment timing.
* Preserved all panel-data functionality from version 0.1.0.

## Implementation notes

* Repeated cross-section ATT(g,t) uses the Wooldridge (2023) pooled QMLE
  with a treatment-by-period interaction. The doubly-robust variant
  augments this with inverse probability weighting on the estimated
  propensity score.
* Sampling weights (when supplied via `weightsname`) are used throughout:
  the outcome regression, the propensity score model, and the pooled
  QMLE. They are multiplied with the IPW factor in the doubly-robust path.
* Analytical SEs for the RCS path use `sandwich::vcovCL` when
  `cluster_var` is supplied, `sandwich::vcovHC` (HC1) otherwise. Panel
  SEs continue to use the influence-function approach from v0.1.0;
  set `boot = TRUE` for fully clustered panel inference.
* The bootstrap automatically resamples whole clusters when
  `cluster_var` is provided, units when `data_type = "panel"`, or
  individual rows when `data_type = "repeated_cross_section"` without
  clustering.
* The compiled C++ helpers from v0.1.0 have been replaced with
  equivalent pure-R implementations, eliminating the Rcpp dependency.
  The package no longer requires compilation.

## Bug fixes

* None.

---

# NonlinearDiD 0.1.0

* Initial CRAN release.
* Panel-data ATT(g,t) estimation under logit, probit, Poisson, negative
  binomial, and linear outcome models.
* Doubly-robust estimator combining outcome regression and propensity
  score weighting.
* `nonlinear_aggte()`: event-study, group, calendar, and simple
  aggregations.
* `nonlinear_pretest()`: joint chi-squared and individual pre-trend
  tests.
* `nonlinear_bounds()`: Manski and parallel-trends-constrained bounds.
* `binary_did_logit()`, `binary_did_probit()`, `binary_did_dr()`: 2x2
  binary DiD estimators.
* `count_did_poisson()`: Poisson QMLE DiD.
* `odds_ratio_did()`: scale-free odds-ratio DiD.
* `sim_binary_panel()`, `sim_count_panel()`: simulation utilities.
