## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  fig.width = 7, fig.height = 4, warning = FALSE, message = FALSE
)
library(NonlinearDiD)
library(ggplot2)
set.seed(101)

## ----demo_scale---------------------------------------------------------------
# Demonstrate scale sensitivity
p_ctrl_pre  <- 0.30; p_ctrl_post <- 0.40
p_treat_pre <- 0.25; p_treat_post <- 0.35

cat("=== Probability Scale ===\n")
cat("Control change:  ", round(p_ctrl_post  - p_ctrl_pre, 4),  "\n")
cat("Treated change:  ", round(p_treat_post - p_treat_pre, 4), "\n")
cat("DiD (prob):      ", round((p_treat_post - p_treat_pre) - (p_ctrl_post - p_ctrl_pre), 4), "\n\n")

cat("=== Log-Odds (Logit) Scale ===\n")
cat("Control change:  ", round(qlogis(p_ctrl_post)  - qlogis(p_ctrl_pre), 4),  "\n")
cat("Treated change:  ", round(qlogis(p_treat_post) - qlogis(p_treat_pre), 4), "\n")
cat("DiD (logit):     ", round((qlogis(p_treat_post) - qlogis(p_treat_pre)) -
                               (qlogis(p_ctrl_post) - qlogis(p_ctrl_pre)), 4), "\n\n")

cat("=== Probit Scale ===\n")
cat("Control change:  ", round(qnorm(p_ctrl_post)  - qnorm(p_ctrl_pre), 4),  "\n")
cat("Treated change:  ", round(qnorm(p_treat_post) - qnorm(p_treat_pre), 4), "\n")
cat("DiD (probit):    ", round((qnorm(p_treat_post) - qnorm(p_treat_pre)) -
                               (qnorm(p_ctrl_post) - qnorm(p_ctrl_pre)), 4), "\n")

## ----severity-----------------------------------------------------------------
# Show severity across baseline probability values
baseline_probs <- seq(0.05, 0.45, by = 0.05)
delta_p        <- 0.10  # same additive change for both groups

severity_df <- do.call(rbind, lapply(baseline_probs, function(p0) {
  p1 <- p0 + delta_p
  # Parallel in prob => same change
  # Logit DiD if treated has different baseline (p0 - 0.05)
  p0_treat <- max(p0 - 0.05, 0.02)
  p1_treat <- p0_treat + delta_p

  logit_did <- (qlogis(p1_treat) - qlogis(p0_treat)) -
               (qlogis(p1)       - qlogis(p0))

  data.frame(
    baseline_ctrl  = p0,
    baseline_treat = p0_treat,
    logit_did      = logit_did
  )
}))

cat("Logit-scale DiD when true probability DiD = 0:\n")
print(severity_df, digits = 3, row.names = FALSE)
cat("\nLarger deviations at low baseline probabilities.\n")

## ----simulation---------------------------------------------------------------
# DGP: parallel trends on logit scale
dat <- sim_binary_panel(
  n            = 1000,
  nperiods     = 8,
  prop_treated = 0.5,
  n_cohorts    = 3,
  true_att     = c(0.20, 0.35, 0.25),
  base_prob    = 0.20,   # low baseline: nonlinearity matters most
  unit_fe_sd   = 0.5,
  seed         = 42
)

cat("Baseline outcome rate (untreated, pre-period):",
    round(mean(dat$y[dat$D == 0 & dat$period == 1]), 3), "\n")
cat("True ATTs (avg):", round(mean(c(0.20, 0.35, 0.25)), 3), "\n\n")

## ----fit_both-----------------------------------------------------------------
# Logit DiD
res_logit <- nonlinear_attgt(
  dat, "y", "period", "id", "g",
  outcome_model = "logit",
  control_group = "nevertreated"
)

# Linear DiD
res_linear <- nonlinear_attgt(
  dat, "y", "period", "id", "g",
  outcome_model = "linear",
  control_group = "nevertreated"
)

# Aggregate
agg_logit  <- nonlinear_aggte(res_logit,  type = "dynamic")
agg_linear <- nonlinear_aggte(res_linear, type = "dynamic")

cat("=== Overall ATT ===\n")
cat("Linear DiD: ", round(agg_linear$overall_att, 4), "\n")
cat("Logit DiD:  ", round(agg_logit$overall_att,  4), "\n")
cat("True ATT:   ", round(mean(c(0.20, 0.35, 0.25)), 4), "\n")

## ----pretrends_both-----------------------------------------------------------
# Test on logit scale
pt_logit <- nonlinear_pretest(res_logit, plot = FALSE)
cat("Pre-trends test (logit scale):\n")
cat("  Joint p-value:", round(pt_logit$joint_pval, 4), "\n\n")

# Test on linear scale
pt_linear <- nonlinear_pretest(res_linear, plot = FALSE)
cat("Pre-trends test (linear scale):\n")
cat("  Joint p-value:", round(pt_linear$joint_pval, 4), "\n\n")

cat("Note: If true DGP is logit-scale parallel trends, the linear-scale\n")
cat("pre-trends test may spuriously reject due to functional form.\n")

