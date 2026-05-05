library(testthat)
library(NonlinearDiD)

# ============================================================
# Test: Data Simulation
# ============================================================

test_that("sim_binary_panel produces valid output", {
  dat <- sim_binary_panel(n = 200, nperiods = 6, seed = 1)
  
  expect_s3_class(dat, "data.frame")
  expect_true(all(c("id", "period", "y", "g", "D") %in% names(dat)))
  expect_equal(nrow(dat), 200 * 6)
  expect_true(all(dat$y %in% c(0, 1)))
  expect_true(all(dat$period %in% 1:6))
  expect_true(all(dat$D %in% c(0, 1)))
})

test_that("sim_binary_panel respects prop_treated", {
  dat <- sim_binary_panel(n = 500, nperiods = 4, prop_treated = 0.4, seed = 10)
  units <- dat[dat$period == 1, ]
  prop <- mean(units$g > 0)
  # Allow tolerance: n_cohorts * n_per_cohort might not be exact
  expect_true(abs(prop - 0.4) < 0.15)
})

test_that("sim_count_panel produces valid count output", {
  dat <- sim_count_panel(n = 100, nperiods = 4, seed = 2)
  expect_s3_class(dat, "data.frame")
  expect_true(all(dat$y >= 0))
  expect_true(is.integer(dat$y) || all(dat$y == floor(dat$y)))
})

# ============================================================
# Test: 2x2 Binary DiD Estimators
# ============================================================

test_that("binary_did_logit returns valid estimates", {
  dat <- sim_binary_panel(n = 300, nperiods = 4, seed = 3)
  dat2 <- dat[dat$period %in% c(2, 3), ]
  
  res <- binary_did_logit(dat2, "y", "period", "id",
                           treat_period = 3, control_period = 2, gname = "g")
  
  expect_s3_class(res, "binary_did_logit")
  expect_true(is.numeric(res$att_link))
  expect_false(is.na(res$att_link))
  expect_true(is.numeric(res$att_ape))
  expect_true(is.numeric(res$se_link))
  expect_true(res$se_link > 0)
  expect_true(res$pval >= 0 && res$pval <= 1)
})

test_that("binary_did_probit returns valid estimates", {
  dat <- sim_binary_panel(n = 300, nperiods = 4, seed = 4)
  dat2 <- dat[dat$period %in% c(2, 4), ]
  
  res <- binary_did_probit(dat2, "y", "period", "id",
                            treat_period = 4, control_period = 2, gname = "g")
  
  expect_s3_class(res, "binary_did_probit")
  expect_false(is.na(res$att_link))
  expect_true(res$se_link > 0)
})

test_that("binary_did_dr returns valid estimate", {
  dat <- sim_binary_panel(n = 400, nperiods = 4, seed = 5)
  dat2 <- dat[dat$period %in% c(2, 4), ]
  
  res <- binary_did_dr(dat2, "y", "period", "id",
                        treat_period = 4, control_period = 2,
                        gname = "g", outcome_model = "logit")
  
  expect_s3_class(res, "binary_did_dr")
  expect_true(is.numeric(res$att))
  expect_false(is.na(res$att))
  expect_true(res$se > 0)
})

test_that("odds_ratio_did returns OR >= 0", {
  dat <- sim_binary_panel(n = 400, nperiods = 4, seed = 6)
  dat2 <- dat[dat$period %in% c(2, 4), ]
  
  res <- odds_ratio_did(dat2, "y", "period", "id",
                         treat_period = 4, control_period = 2, gname = "g")
  
  expect_s3_class(res, "odds_ratio_did")
  expect_true(res$or_did > 0)
  expect_true(is.numeric(res$log_or_did))
  expect_equal(exp(res$log_or_did), res$or_did, tolerance = 1e-8)
})

# ============================================================
# Test: Count DiD
# ============================================================

test_that("count_did_poisson returns valid rate ratio", {
  dat <- sim_count_panel(n = 300, nperiods = 6, seed = 7)
  dat2 <- dat[dat$period %in% c(2, 4), ]
  
  res <- count_did_poisson(dat2, "y", "period", "id",
                            treat_period = 4, control_period = 2, gname = "g")
  
  expect_s3_class(res, "count_did_poisson")
  expect_true(res$rate_ratio > 0)
  expect_equal(exp(res$att_log_rr), res$rate_ratio, tolerance = 1e-8)
  expect_true(res$ci_lo_rr < res$rate_ratio)
  expect_true(res$ci_hi_rr > res$rate_ratio)
})

# ============================================================
# Test: ATT(g,t) Estimation
# ============================================================

test_that("nonlinear_attgt produces a valid nonlinear_attgt object", {
  dat <- sim_binary_panel(n = 300, nperiods = 6, n_cohorts = 2, seed = 8)
  
  res <- nonlinear_attgt(
    data          = dat,
    yname         = "y",
    tname         = "period",
    idname        = "id",
    gname         = "g",
    outcome_model = "logit",
    control_group = "nevertreated"
  )
  
  expect_s3_class(res, "nonlinear_attgt")
  expect_true(is.data.frame(res$attgt))
  expect_true(all(c("group", "time", "att", "post") %in% names(res$attgt)))
  expect_true(nrow(res$attgt) > 0)
})

test_that("nonlinear_attgt errors on missing columns", {
  dat <- sim_binary_panel(n = 100, nperiods = 4, seed = 9)
  expect_error(
    nonlinear_attgt(dat, "y", "period", "id", "BADCOL"),
    "Missing columns"
  )
})

test_that("nonlinear_attgt errors with no treated units", {
  dat <- sim_binary_panel(n = 100, nperiods = 4, seed = 10)
  dat$g <- 0L  # remove all treatment
  expect_error(
    nonlinear_attgt(dat, "y", "period", "id", "g"),
    "No treated units"
  )
})

test_that("nonlinear_attgt works with probit model", {
  dat <- sim_binary_panel(n = 300, nperiods = 5, n_cohorts = 2, seed = 11)
  
  res <- nonlinear_attgt(dat, "y", "period", "id", "g",
                          outcome_model = "probit")
  expect_s3_class(res, "nonlinear_attgt")
})

test_that("nonlinear_attgt works with poisson on count data", {
  dat <- sim_count_panel(n = 200, nperiods = 5, n_cohorts = 2, seed = 12)
  res <- nonlinear_attgt(dat, "y", "period", "id", "g",
                          outcome_model = "poisson")
  expect_s3_class(res, "nonlinear_attgt")
})

# ============================================================
# Test: Aggregation
# ============================================================

test_that("nonlinear_aggte produces valid dynamic aggregation", {
  dat <- sim_binary_panel(n = 300, nperiods = 6, n_cohorts = 2, seed = 13)
  res <- nonlinear_attgt(dat, "y", "period", "id", "g",
                          outcome_model = "logit")
  
  agg <- nonlinear_aggte(res, type = "dynamic")
  
  expect_s3_class(agg, "nonlinear_aggte")
  expect_equal(agg$type, "dynamic")
  expect_true(is.numeric(agg$overall_att))
  expect_true(is.data.frame(agg$agg))
  expect_true("label" %in% names(agg$agg))
})

test_that("nonlinear_aggte works for all types", {
  dat <- sim_binary_panel(n = 300, nperiods = 6, n_cohorts = 2, seed = 14)
  res <- nonlinear_attgt(dat, "y", "period", "id", "g",
                          outcome_model = "logit")
  
  for (tp in c("dynamic", "group", "calendar", "simple")) {
    agg <- nonlinear_aggte(res, type = tp)
    expect_s3_class(agg, "nonlinear_aggte")
    expect_equal(agg$type, tp)
  }
})

# ============================================================
# Test: Bounds
# ============================================================

test_that("nonlinear_bounds returns lb <= ub", {
  dat <- sim_binary_panel(n = 200, nperiods = 5, n_cohorts = 2, seed = 15)
  bounds <- nonlinear_bounds(dat, "y", "period", "id", "g",
                              bound_type = "manski")
  
  valid <- bounds[!is.na(bounds$lb) & !is.na(bounds$ub), ]
  expect_true(all(valid$lb <= valid$ub))
})

test_that("nonlinear_bounds PT identification yields equal lb and ub", {
  dat <- sim_binary_panel(n = 300, nperiods = 5, n_cohorts = 2, seed = 16)
  bounds <- nonlinear_bounds(dat, "y", "period", "id", "g",
                              bound_type = "pt_only")
  
  valid <- bounds[!is.na(bounds$lb) & !is.na(bounds$ub), ]
  # Under PT alone, point identified => lb == ub
  expect_true(all(abs(valid$lb - valid$ub) < 1e-10))
})

# ============================================================
# Test: Pre-trends Test
# ============================================================

test_that("nonlinear_pretest returns valid test object", {
  dat <- sim_binary_panel(n = 400, nperiods = 8, n_cohorts = 2, seed = 17)
  res <- nonlinear_attgt(dat, "y", "period", "id", "g",
                          outcome_model = "logit")
  
  pt  <- nonlinear_pretest(res, plot = FALSE)
  
  expect_s3_class(pt, "nonlinear_pretest")
  expect_true(is.numeric(pt$joint_stat))
  expect_true(is.numeric(pt$joint_pval))
  expect_true(pt$joint_pval >= 0 && pt$joint_pval <= 1)
  expect_true(is.character(pt$conclusion))
})

# ============================================================
# Test: Print/Summary (no errors)
# ============================================================

test_that("print and summary methods run without errors", {
  dat <- sim_binary_panel(n = 200, nperiods = 5, n_cohorts = 2, seed = 18)
  res <- nonlinear_attgt(dat, "y", "period", "id", "g",
                          outcome_model = "logit")
  
  expect_output(print(res))
  expect_output(summary(res))
  
  agg <- nonlinear_aggte(res, type = "dynamic")
  expect_output(print(agg))
  expect_output(summary(agg))
})

test_that("print methods for 2x2 estimators run without errors", {
  dat <- sim_binary_panel(n = 200, nperiods = 4, seed = 19)
  dat2 <- dat[dat$period %in% c(2, 3), ]
  
  expect_output(print(binary_did_logit(dat2, "y", "period", "id", 3, 2, gname = "g")))
  expect_output(print(binary_did_probit(dat2, "y", "period", "id", 3, 2, gname = "g")))
  expect_output(print(binary_did_dr(dat2, "y", "period", "id", 3, 2, gname = "g")))
  expect_output(print(odds_ratio_did(dat2, "y", "period", "id", 3, 2, gname = "g")))
})
