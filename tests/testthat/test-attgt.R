test_that("v0.1.0 panel syntax still works (logit, default DR)", {
  set.seed(1)
  dat <- sim_binary_panel(n = 300, nperiods = 4, prop_treated = 0.5)
  res <- nonlinear_attgt(
    data          = dat,
    yname         = "y",
    tname         = "period",
    idname        = "id",
    gname         = "g",
    outcome_model = "logit"
  )
  expect_s3_class(res, "nonlinear_attgt")
  expect_true("att" %in% names(res$attgt))
  expect_equal(res$args$data_type, "panel")
})

test_that("idname is required for panel data", {
  set.seed(2)
  dat <- sim_binary_panel(n = 200, nperiods = 4, prop_treated = 0.5)
  expect_error(
    nonlinear_attgt(dat, yname = "y", tname = "period",
                    gname = "g", outcome_model = "logit"),
    "idname"
  )
})

test_that("repeated cross-section runs without idname", {
  set.seed(3)
  rcs <- sim_binary_rcs(n_per_period = 250, nperiods = 5, prop_treated = 0.5)
  res <- nonlinear_attgt(
    data          = rcs,
    yname         = "y",
    tname         = "period",
    gname         = "g",
    outcome_model = "logit",
    data_type     = "repeated_cross_section"
  )
  expect_s3_class(res, "nonlinear_attgt")
  expect_equal(res$args$data_type, "repeated_cross_section")
  # At least some post-treatment ATTs should be estimable
  expect_true(any(!is.na(res$attgt$att[res$attgt$post])))
})

test_that("weightsname is accepted and uses the weight column", {
  set.seed(4)
  rcs <- sim_binary_rcs(n_per_period = 250, nperiods = 5, prop_treated = 0.5)
  rcs$w <- runif(nrow(rcs), 0.5, 1.5)
  res <- nonlinear_attgt(
    data          = rcs,
    yname         = "y",
    tname         = "period",
    gname         = "g",
    outcome_model = "logit",
    data_type     = "repeated_cross_section",
    weightsname   = "w"
  )
  expect_s3_class(res, "nonlinear_attgt")
  expect_equal(res$args$weightsname, "w")
})

test_that("cluster_var produces non-NA standard errors via vcovCL", {
  set.seed(5)
  rcs <- sim_binary_rcs(n_per_period = 300, nperiods = 5, prop_treated = 0.5)
  # Synthetic cluster: 10 groups
  rcs$cl <- sample(1:10, nrow(rcs), replace = TRUE)
  res <- nonlinear_attgt(
    data          = rcs,
    yname         = "y",
    tname         = "period",
    gname         = "g",
    outcome_model = "logit",
    data_type     = "repeated_cross_section",
    cluster_var   = "cl"
  )
  expect_s3_class(res, "nonlinear_attgt")
  expect_true(any(!is.na(res$attgt$se)))
})

test_that("APE estimand returns probability-scale effects", {
  set.seed(6)
  rcs <- sim_binary_rcs(n_per_period = 250, nperiods = 5, prop_treated = 0.5)
  res <- nonlinear_attgt(
    data          = rcs,
    yname         = "y",
    tname         = "period",
    gname         = "g",
    outcome_model = "logit",
    data_type     = "repeated_cross_section",
    estimand      = "ape"
  )
  expect_equal(res$args$estimand, "ape")
  # APE estimates should be on probability scale: |ape| <= 1
  expect_true(all(abs(res$attgt$att[!is.na(res$attgt$att)]) <= 1))
})

test_that("binary_did_dr 2x2 estimator still works (regression test)", {
  set.seed(7)
  dat <- sim_binary_panel(n = 300, nperiods = 4, prop_treated = 0.5)
  dat2 <- dat[dat$period %in% c(2, 3), ]
  res <- binary_did_dr(dat2, "y", "period", "id", 3, 2,
                       gname = "g", outcome_model = "logit")
  expect_s3_class(res, "binary_did_dr")
  expect_true(is.numeric(res$att))
})

test_that("doubly_robust = FALSE works for RCS", {
  set.seed(8)
  rcs <- sim_binary_rcs(n_per_period = 250, nperiods = 5, prop_treated = 0.5)
  res <- nonlinear_attgt(
    data          = rcs,
    yname         = "y",
    tname         = "period",
    gname         = "g",
    outcome_model = "logit",
    data_type     = "repeated_cross_section",
    doubly_robust = FALSE
  )
  expect_s3_class(res, "nonlinear_attgt")
})
