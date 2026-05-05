## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width  = 7,
  fig.height = 4.5,
  warning    = FALSE,
  message    = FALSE
)
library(NonlinearDiD)
library(ggplot2)
set.seed(42)

## ----simulate-----------------------------------------------------------------
dat <- sim_binary_panel(
  n            = 800,
  nperiods     = 8,
  prop_treated = 0.6,
  n_cohorts    = 3,
  true_att     = c(0.15, 0.25, 0.20),  # heterogeneous treatment effects
  base_prob    = 0.3,
  seed         = 42
)

head(dat)
cat("Treatment cohorts:", table(dat$g[dat$period == 1]), "\n")
cat("Baseline outcome rate:", round(mean(dat$y[dat$D == 0]), 3), "\n")

## ----attgt--------------------------------------------------------------------
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

summary(res)

## ----eventstudy, fig.cap="Event-study plot for binary outcome DiD"------------
agg_dynamic <- nonlinear_aggte(res, type = "dynamic")
print(agg_dynamic)
plot(agg_dynamic)

## ----groupatt-----------------------------------------------------------------
agg_group <- nonlinear_aggte(res, type = "group")
print(agg_group)

## ----pretest------------------------------------------------------------------
pt <- nonlinear_pretest(res, plot = FALSE)
print(pt)

## ----ordata-------------------------------------------------------------------
# Use simple 2-period case for clarity
dat2 <- dat[dat$period %in% c(3, 5), ]

res_or <- odds_ratio_did(
  data           = dat2,
  yname          = "y",
  tname          = "period",
  idname         = "id",
  treat_period   = 5,
  control_period = 3,
  gname          = "g"
)

print(res_or)

## ----countdata----------------------------------------------------------------
count_dat <- sim_count_panel(
  n            = 600,
  nperiods     = 6,
  prop_treated = 0.5,
  true_rr      = c(1.5, 2.0, 1.3),  # rate ratios per cohort
  base_rate    = 10,
  seed         = 7
)

summary(count_dat$y)

## ----countdid-----------------------------------------------------------------
# Staggered Poisson DiD
res_count <- nonlinear_attgt(
  data          = count_dat,
  yname         = "y",
  tname         = "period",
  idname        = "id",
  gname         = "g",
  outcome_model = "poisson",
  estimand      = "att",
  control_group = "nevertreated"
)

agg_count <- nonlinear_aggte(res_count, type = "dynamic")
plot(agg_count)

## ----pois2x2------------------------------------------------------------------
# Simple 2x2 Poisson DiD (Wooldridge QMLE)
count_sub <- count_dat[count_dat$period %in% c(2, 4), ]
res_pois  <- count_did_poisson(
  count_sub,
  yname          = "y",
  tname          = "period",
  idname         = "id",
  treat_period   = 4,
  control_period = 2,
  gname          = "g"
)
print(res_pois)

## ----drlogit------------------------------------------------------------------
res_dr <- binary_did_dr(
  data           = dat[dat$period %in% c(3, 5), ],
  yname          = "y",
  tname          = "period",
  idname         = "id",
  treat_period   = 5,
  control_period = 3,
  gname          = "g",
  xformla        = ~ x1 + x2,
  outcome_model  = "logit"
)
print(res_dr)

## ----bounds-------------------------------------------------------------------
bounds <- nonlinear_bounds(
  data          = dat,
  yname         = "y",
  tname         = "period",
  idname        = "id",
  gname         = "g",
  bound_type    = "manski",  # widest (no assumptions)
  control_group = "nevertreated"
)

# Show post-treatment bounds
post_bounds <- bounds[bounds$post, ]
head(post_bounds[, c("group", "time", "lb", "ub", "identified")])

## ----comparison, eval = FALSE-------------------------------------------------
# # Linear comparison (for illustration)
# res_linear <- nonlinear_attgt(
#   data          = dat,
#   yname         = "y",
#   tname         = "period",
#   idname        = "id",
#   gname         = "g",
#   outcome_model = "linear",
#   estimand      = "att"
# )
# 
# agg_lin  <- nonlinear_aggte(res_linear, type = "dynamic")
# agg_logit <- nonlinear_aggte(res, type = "dynamic")
# 
# # The two produce different estimates when baseline rates are moderate
# cat("Linear overall ATT:",  round(agg_lin$overall_att, 4), "\n")
# cat("Logit overall ATT:",   round(agg_logit$overall_att, 4), "\n")
# cat("True ATT (avg):    ",  round(mean(c(0.15, 0.25, 0.20)), 4), "\n")

