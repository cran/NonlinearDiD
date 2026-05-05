# NonlinearDiD 0.1.0

## Initial Release

### New Features

* `nonlinear_attgt()`: Core estimator for group-time ATT(g,t) under logit,
  probit, Poisson, negative binomial, and linear outcome models with staggered
  treatment adoption.

* `nonlinear_aggte()`: Aggregation into event-study (dynamic), group-level,
  calendar-time, and overall ATT estimates.

* `nonlinear_pretest()`: Pre-treatment parallel trends test with joint
  chi-squared test and HonestDiD-style sensitivity analysis.

* `binary_did_logit()` / `binary_did_probit()`: Simple 2×2 DiD with binary
  outcomes on the log-odds / probit scale, with APE reporting.

* `binary_did_dr()`: Doubly-robust binary DiD combining logit outcome
  regression with inverse probability weighting.

* `count_did_poisson()`: Poisson QMLE DiD for count outcomes following
  Wooldridge (2023), reporting rate ratios.

* `odds_ratio_did()`: Odds-ratio DiD estimator (scale-free, symmetric).

* `nonlinear_bounds()`: Nonparametric Manski bounds and PT-restricted bounds
  for binary outcomes.

* `sim_binary_panel()` / `sim_count_panel()`: Data-generating processes for
  simulation studies with staggered treatment and heterogeneous effects.

* S3 methods: `print()`, `summary()`, `plot()` for all main object classes.

* Rcpp-accelerated bootstrap weight generation and DR score computation.

### Methodological Notes

This is version 0.1.0 — an initial implementation of a methodology that is
actively being developed in the econometrics literature. The core identification
arguments follow Roth & Sant'Anna (2023) and Wooldridge (2023). Standard errors
are based on influence function / sandwich estimators.

Known limitations:
- Simultaneous confidence bands use a conservative normal approximation;
  exact bands require the multiplier bootstrap (`boot = TRUE`).
- Negative binomial staggered DiD uses approximation; full MLE version
  is planned for v0.2.0.
