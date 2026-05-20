# NonlinearDiD

**Staggered Difference-in-Differences with Nonlinear Outcomes — panel and repeated cross-section data**

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen)](https://github.com/causalfragility-lab/NonlinearDiD)
[![CRAN status](https://www.r-pkg.org/badges/version/NonlinearDiD)](https://CRAN.R-project.org/package=NonlinearDiD)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## The Problem

The Callaway & Sant'Anna (2021) framework for staggered DiD is hugely influential — but it assumes continuous outcomes with linear parallel trends.

For binary outcomes (employed/not, hospitalized/not, defaulted/not), this creates fundamental problems:

| Problem | Why it matters |
|---|---|
| **Scale sensitivity** | Parallel trends in P(Y=1) ≠ parallel trends in log-odds. Pre-trends can appear flat or steep depending on the scale. |
| **Jensen's inequality** | Treatment effects on the probability scale mix the "real" effect with curvature of the CDF. |
| **Heterogeneous baseline rates** | Units with different baseline probabilities will show "spurious" violations of parallel trends even under no treatment effect. |

`NonlinearDiD` extends the CS2021 framework to properly handle logit, probit, Poisson, and negative binomial outcome models — for both **panel** and **repeated cross-section** data.

## What's New in 0.2.0

- **Repeated cross-section support** — different individuals each period (e.g. BRFSS, NHIS, CPS supplements). Set `data_type = "repeated_cross_section"`; `idname` becomes optional.
- **Sampling weights** — pass `weightsname = "wt"` and the weight is threaded through the outcome regression, propensity score, and pooled QMLE.
- **Clustered inference** — pass `cluster_var = "state"` for `sandwich::vcovCL()` analytical SEs and cluster-resampling bootstrap.
- **No more compilation** — the Rcpp helpers from 0.1.0 are now pure R, so installation is one step on every platform.
- **All v0.1.0 functions and arguments preserved** — existing scripts using named arguments continue to work unchanged.

## Installation

```r
# CRAN
install.packages("NonlinearDiD")

# GitHub (development version)
remotes::install_github("causalfragility-lab/NonlinearDiD")
```

## Quick Start: Panel Data

```r
library(NonlinearDiD)

# 1. Simulate staggered binary panel data
dat <- sim_binary_panel(n = 500, nperiods = 8, n_cohorts = 3,
                        prop_treated = 0.5, true_att = 0.25, seed = 42)

# 2. Estimate ATT(g,t) with logistic outcome model
res <- nonlinear_attgt(
  data          = dat,
  yname         = "y",
  tname         = "period",
  idname        = "id",
  gname         = "g",
  xformla       = ~ x1 + x2,
  outcome_model = "logit",
  estimand      = "att",
  control_group = "nevertreated",
  doubly_robust = TRUE
)

# 3. Aggregate into event-study
agg <- nonlinear_aggte(res, type = "dynamic")
plot(agg)

# 4. Pre-treatment parallel trends test
nonlinear_pretest(res)
```

## Quick Start: Repeated Cross-Section Data

```r
library(NonlinearDiD)

# 1. Simulate repeated cross-section binary data
rcs <- sim_binary_rcs(n_per_period = 500, nperiods = 8,
                      prop_treated = 0.5, true_att = 0.3, seed = 7)

# 2. Estimate ATT(g,t) — note the new data_type argument
res <- nonlinear_attgt(
  data          = rcs,
  yname         = "y",
  tname         = "period",
  gname         = "g",
  outcome_model = "logit",
  estimand      = "ape",
  data_type     = "repeated_cross_section",
  control_group = "notyetreated"
)

plot(nonlinear_aggte(res, type = "dynamic"))
```

## Survey-Weighted Real-World Example

Repeated cross-section survey data (e.g. CPS Food Security Supplement) with sampling weights and state-level clustering:

```r
res <- nonlinear_attgt(
  data          = snap_data,
  yname         = "food_insecure",
  tname         = "year",
  gname         = "policy_end_year",
  idname        = "household_id",    # optional — used only as a record ID
  data_type     = "repeated_cross_section",
  outcome_model = "logit",
  estimand      = "ape",
  weightsname   = "survey_weight",
  cluster_var   = "state",
  control_group = "notyetreated"
)

summary(res)
nonlinear_aggte(res, type = "dynamic")
```

## Key Functions

| Function | Description |
|---|---|
| `nonlinear_attgt()` | Main engine: estimates ATT(g,t) for all cohort × time cells. Panel and repeated cross-section. |
| `nonlinear_aggte()` | Aggregates ATT(g,t) into event-study, group, calendar, or overall ATT |
| `nonlinear_pretest()` | Tests pre-treatment parallel trends (joint + individual + HonestDiD) |
| `binary_did_logit()` | Simple 2×2 DiD with logistic outcome |
| `binary_did_probit()` | Simple 2×2 DiD with probit outcome |
| `binary_did_dr()` | Doubly-robust binary DiD |
| `count_did_poisson()` | Poisson QMLE DiD for count outcomes (Wooldridge 2023) |
| `odds_ratio_did()` | Odds-ratio DiD estimator |
| `nonlinear_bounds()` | Nonparametric Manski / PT bounds |
| `sim_binary_panel()` | Simulate binary panel data for testing |
| `sim_binary_rcs()` | Simulate binary repeated cross-section data **(new in 0.2.0)** |
| `sim_count_panel()` | Simulate count panel data for testing |

## Estimands

| Estimand | Scale | When to use |
|---|---|---|
| `"att"` | Link scale (log-odds / probit index / log-count) | Compare with linear DiD coefficient |
| `"ape"` | Probability scale | What practitioners usually report |
| `"odds_ratio"` | Multiplicative | Scale-free; natural for 2×2 tables |

## Panel vs Repeated Cross-Section

`NonlinearDiD` supports staggered difference-in-differences designs with nonlinear outcomes for both panel and repeated cross-section data.

- **Panel data** (`data_type = "panel"`, default): units are followed over time and `idname` identifies repeated observations. Estimation uses within-unit outcome changes following Callaway & Sant'Anna (2021).

- **Repeated cross-section data** (`data_type = "repeated_cross_section"`): observations are independent within each time period. `idname` is optional and may identify survey records or households, but the estimator does not require the same units to appear across periods. Estimation uses pooled quasi-maximum likelihood approaches motivated by Wooldridge (2023), with an optional IPW-augmented doubly-robust variant.

## Outcome Models

| `outcome_model` | Parallel Trends Assumption | Outcome Type |
|---|---|---|
| `"logit"` | Parallel in log-odds | Binary (0/1) |
| `"probit"` | Parallel in probit index | Binary (0/1) |
| `"poisson"` | Parallel in log-count | Count (≥ 0) |
| `"negbin"` | Parallel in log-count | Overdispersed count |
| `"linear"` | Parallel in mean (LPM) | Continuous / binary |

## References

- Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200–230. <https://doi.org/10.1016/j.jeconom.2020.12.001>
- Roth, J., & Sant'Anna, P. H. C. (2023). When is parallel trends sensitive to functional form? *Econometrica*, 91(2), 737–747. <https://doi.org/10.3982/ECTA19402>
- Wooldridge, J. M. (2023). Simple approaches to nonlinear difference-in-differences with panel data. *The Econometrics Journal*, 26(3), C31–C66. <https://doi.org/10.1093/ectj/utad016>
- Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics*, 219(1), 101–122. <https://doi.org/10.1016/j.jeconom.2020.06.003>
- Manski, C. F. (1990). Nonparametric bounds on treatment effects. *American Economic Review*, 80(2), 319–323.

## Contributing

This package addresses an active research frontier. Contributions, bug reports, and methodological suggestions are welcome — please open an issue or pull request on GitHub.

## License

MIT © 2026 Subir Hait
