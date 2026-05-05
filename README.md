# NonlinearDiD <img src="man/figures/logo.png" align="right" height="139" />

> **Staggered Difference-in-Differences with Nonlinear Outcomes**

[![R-CMD-check](https://github.com/causalfragility-lab/NonlinearDiD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/causalfragility-lab/NonlinearDiD/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/NonlinearDiD)](https://CRAN.R-project.org/package=NonlinearDiD)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## The Problem

The **Callaway & Sant'Anna (2021)** framework for staggered DiD is hugely
influential — but it assumes **continuous outcomes with linear parallel trends**.

For binary outcomes (employed/not, hospitalized/not, defaulted/not), this
creates fundamental problems:

| Problem | Why it matters |
|---------|----------------|
| **Scale sensitivity** | Parallel trends in P(Y=1) ≠ parallel trends in log-odds. Pre-trends can appear flat or steep depending on the scale. |
| **Jensen's inequality** | Treatment effects on probability scale mix the "real" effect with curvature of the CDF. |
| **Heterogeneous baseline rates** | Units with different baseline probabilities will show "spurious" violations of parallel trends even under no treatment effect. |

**NonlinearDiD** solves all three by extending the CS2021 framework to
properly handle logit, probit, Poisson, and negative binomial outcome models.

---

## Installation

```r
# CRAN
install.packages("NonlinearDiD")

# GitHub (development version)
remotes::install_github("causalfragility-lab/NonlinearDiD")
```

---

## Quick Start

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
  xformla       = ~ x1 + x2,      # covariates
  outcome_model = "logit",         # logit / probit / poisson / negbin
  estimand      = "att",           # att / odds_ratio / ape
  control_group = "nevertreated",  # nevertreated / notyetreated
  doubly_robust = TRUE             # doubly-robust estimator
)

# 3. Aggregate into event-study
agg <- nonlinear_aggte(res, type = "dynamic")
plot(agg)

# 4. Pre-treatment parallel trends test
nonlinear_pretest(res)
```

---

## Key Functions

| Function | Description |
|----------|-------------|
| `nonlinear_attgt()` | Main engine: estimates ATT(g,t) for all cohort × time cells |
| `nonlinear_aggte()` | Aggregates ATT(g,t) into event-study, group, calendar, or overall ATT |
| `nonlinear_pretest()` | Tests pre-treatment parallel trends (joint + individual + HonestDiD) |
| `binary_did_logit()` | Simple 2×2 DiD with logistic outcome |
| `binary_did_probit()` | Simple 2×2 DiD with probit outcome |
| `binary_did_dr()` | Doubly-robust binary DiD |
| `count_did_poisson()` | Poisson QMLE DiD for count outcomes (Wooldridge 2023) |
| `odds_ratio_did()` | Odds-ratio DiD estimator |
| `nonlinear_bounds()` | Nonparametric Manski/PT bounds |
| `sim_binary_panel()` | Simulate binary panel data for testing |
| `sim_count_panel()` | Simulate count panel data for testing |

---

## Estimands

- **ATT** (`estimand = "att"`): Average treatment effect on the treated, on
  the linear probability / log-odds / log-count scale depending on `outcome_model`.
- **APE** (`estimand = "ape"`): Average partial effect on the probability
  scale (for logit/probit) — what practitioners usually want.
- **Odds-ratio** (`estimand = "odds_ratio"`): Multiplicative odds ratio DiD;
  invariant to group labeling; natural for 2×2 tables.

---

## Methods Reference

### Outcome Models

| `outcome_model` | Parallel Trends Assumption | Outcome Type |
|-----------------|---------------------------|--------------|
| `"logit"` | Parallel in log-odds | Binary (0/1) |
| `"probit"` | Parallel in probit index | Binary (0/1) |
| `"poisson"` | Parallel in log-count | Count (≥ 0) |
| `"negbin"` | Parallel in log-count | Overdispersed count |
| `"linear"` | Parallel in mean (LPM) | Continuous / binary |

---

## References

- Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with
  multiple time periods. *Journal of Econometrics*, 225(2), 200–230.

- Roth, J., & Sant'Anna, P. H. C. (2023). When is parallel trends sensitive to
  functional form? *Econometrica*, 91(2), 737–747.

- Wooldridge, J. M. (2023). Simple approaches to nonlinear
  difference-in-differences with panel data. *The Econometrics Journal*, 26(3), 31–66.

- Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust difference-in-differences
  estimators. *Journal of Econometrics*, 219(1), 101–122.

- Manski, C. F. (1990). Nonparametric bounds on treatment effects.
  *American Economic Review*, 80(2), 319–323.

---

## Contributing

This package addresses an active research frontier. Contributions,
bug reports, and methodological suggestions are welcome — please open
an issue or pull request on
[GitHub](https://github.com/causalfragility-lab/NonlinearDiD/issues).

---

## License

MIT © 2026 Subir Hait
