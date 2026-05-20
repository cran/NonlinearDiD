#' @title Nonlinear Staggered DiD: Group-Time ATT Estimation
#'
#' @description
#' Computes group-time average treatment effects on the treated (ATT(g,t)) for
#' staggered difference-in-differences designs with nonlinear outcomes.
#' Supports both \strong{panel} data (same units across periods) and
#' \strong{repeated cross-section} (RCS) data (independent samples per period).
#'
#' For panel data the package follows Callaway & Sant'Anna (2021) and uses
#' within-unit outcome changes to estimate counterfactual trends. For repeated
#' cross-sections it uses the Wooldridge (2023) pooled QMLE with a
#' treatment-by-period interaction (non-DR) or an IPW-augmented version
#' (doubly-robust). Both modes optionally accept sampling weights and a
#' clustering variable.
#'
#' @param data A data frame in long format.
#' @param yname Character. Outcome variable column.
#' @param tname Character. Time period column.
#' @param gname Character. Treatment cohort column (the period when a
#'   unit/group first receives treatment; 0 or Inf for never-treated).
#' @param idname Character or \code{NULL}. Unit identifier column.
#'   Required for \code{data_type = "panel"}. Optional for
#'   \code{data_type = "repeated_cross_section"}.
#' @param data_type Character. \code{"panel"} (default) or
#'   \code{"repeated_cross_section"}.
#' @param weightsname Character or \code{NULL}. Column name of sampling
#'   weights (e.g. survey design weights). Used in all model fits
#'   (outcome regression, propensity score, pooled QMLE) when supplied.
#'   Default \code{NULL} (equal weights).
#' @param cluster_var Character or \code{NULL}. Column name to cluster
#'   standard errors on (e.g. \code{"state"}). Analytical SEs use
#'   \code{sandwich::vcovCL()} and the bootstrap resamples whole
#'   clusters. Default \code{NULL} (HC1 robust SEs / row resampling).
#' @param xformla A one-sided formula for covariates (e.g. \code{~ x1 + x2}).
#'   Default \code{~ 1}.
#' @param outcome_model Character. One of \code{"logit"}, \code{"probit"},
#'   \code{"poisson"}, \code{"negbin"}, \code{"linear"}.
#' @param estimand Character. \code{"att"} (default), \code{"ape"} (average
#'   partial effect on probability scale), or \code{"odds_ratio"}.
#' @param control_group Character. \code{"nevertreated"} (default) or
#'   \code{"notyetreated"}.
#' @param doubly_robust Logical. Use the doubly-robust estimator. Default TRUE.
#' @param boot Logical. Bootstrap inference. Default FALSE.
#' @param nboot Integer. Bootstrap iterations. Default 999.
#' @param boot_type Character. \code{"multiplier"} or \code{"empirical"}.
#' @param alpha Numeric. Significance level. Default 0.05.
#' @param parallel Logical. Parallel bootstrap. Default FALSE.
#' @param pl_cores Integer. Cores for parallel bootstrap.
#' @param anticipation Integer. Periods of anticipation allowed. Default 0.
#'
#' @return An object of class \code{nonlinear_attgt}.
#'
#' @references
#' Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with
#' multiple time periods. \emph{Journal of Econometrics}, 225(2), 200-230.
#'
#' Wooldridge, J. M. (2023). Simple approaches to nonlinear
#' difference-in-differences with panel data. \emph{The Econometrics Journal}, 26(3).
#'
#' Roth, J., & Sant'Anna, P. H. C. (2023). When is parallel trends sensitive
#' to functional form? \emph{Econometrica}, 91(2), 737-747.
#'
#' Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust
#' difference-in-differences estimators. \emph{Journal of Econometrics},
#' 219(1), 101-122.
#'
#' @examples
#' # ---- Panel example (v0.1.0 syntax — unchanged) ----
#' set.seed(42)
#' dat <- sim_binary_panel(n = 500, nperiods = 6, prop_treated = 0.4)
#' result <- nonlinear_attgt(
#'   data = dat, yname = "y", tname = "period",
#'   idname = "id", gname = "g",
#'   outcome_model = "logit"
#' )
#' summary(result)
#'
#' # ---- Repeated cross-section example ----
#' set.seed(7)
#' rcs <- sim_binary_rcs(n_per_period = 400, nperiods = 6, prop_treated = 0.4)
#' res_rcs <- nonlinear_attgt(
#'   data = rcs, yname = "y", tname = "period", gname = "g",
#'   outcome_model = "logit",
#'   data_type = "repeated_cross_section"
#' )
#' summary(res_rcs)
#'
#' @export
nonlinear_attgt <- function(
    data,
    yname,
    tname,
    gname,
    idname        = NULL,
    data_type     = c("panel", "repeated_cross_section"),
    weightsname   = NULL,
    cluster_var   = NULL,
    xformla       = ~1,
    outcome_model = c("logit", "probit", "poisson", "negbin", "linear"),
    estimand      = c("att", "ape", "odds_ratio"),
    control_group = c("nevertreated", "notyetreated"),
    doubly_robust = TRUE,
    boot          = FALSE,
    nboot         = 999,
    boot_type     = c("multiplier", "empirical"),
    alpha         = 0.05,
    parallel      = FALSE,
    pl_cores      = 2L,
    anticipation  = 0L
) {

  outcome_model <- match.arg(outcome_model)
  estimand      <- match.arg(estimand)
  control_group <- match.arg(control_group)
  boot_type     <- match.arg(boot_type)
  data_type     <- match.arg(data_type)

  if (data_type == "panel" && is.null(idname)) {
    stop(paste0(
      "'idname' is required for data_type = 'panel'.\n",
      "If your data are repeated cross-sections (different units each period),\n",
      "use data_type = 'repeated_cross_section' instead."
    ))
  }

  if (estimand == "odds_ratio" && !outcome_model %in% c("logit", "probit", "linear")) {
    stop("'odds_ratio' estimand requires outcome_model = 'logit', 'probit', or 'linear'.")
  }

  required_cols <- c(yname, tname, gname)
  if (data_type == "panel") required_cols <- c(required_cols, idname)
  if (!is.null(weightsname)) required_cols <- c(required_cols, weightsname)
  if (!is.null(cluster_var)) required_cols <- c(required_cols, cluster_var)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }

  data <- as.data.frame(data)
  data[[gname]] <- ifelse(is.na(data[[gname]]) | is.infinite(data[[gname]]),
                          0, data[[gname]])

  # Internal sampling-weight column .sw — always present, defaults to 1
  data$.sw <- if (!is.null(weightsname)) {
    w <- data[[weightsname]]
    if (any(is.na(w)) || any(w < 0)) {
      stop("'weightsname' must contain non-negative, non-missing values.")
    }
    as.numeric(w)
  } else {
    rep(1.0, nrow(data))
  }

  tlist <- sort(unique(data[[tname]]))
  glist <- sort(setdiff(unique(data[[gname]]), 0))

  if (length(glist) == 0) stop("No treated units found (all gname values are 0 or missing).")
  if (length(tlist) < 2)  stop("Need at least 2 time periods.")

  attgt_list <- vector("list", length(glist) * length(tlist))
  idx <- 0L

  for (g in glist) {
    for (t in tlist) {
      idx <- idx + 1L
      pre_period <- g - 1L - anticipation

      attgt_list[[idx]] <- .compute_attgt_single(
        data          = data,
        yname         = yname,
        tname         = tname,
        idname        = idname,
        gname         = gname,
        xformla       = xformla,
        g             = g,
        t             = t,
        base_period   = pre_period,
        outcome_model = outcome_model,
        estimand      = estimand,
        control_group = control_group,
        doubly_robust = doubly_robust,
        data_type     = data_type
      )
    }
  }

  attgt_df <- do.call(rbind, lapply(attgt_list, function(x) {
    if (is.null(x)) return(NULL)
    data.frame(
      group     = x$g,
      time      = x$t,
      att       = x$att,
      se        = NA_real_,
      post      = x$t >= x$g,
      estimand  = estimand,
      converged = x$converged,
      n_treated = x$n_treated,
      n_control = x$n_control
    )
  }))

  if (boot) {
    attgt_df <- .bootstrap_inference(
      attgt_df = attgt_df, data = data, yname = yname, tname = tname,
      idname = idname, gname = gname, xformla = xformla,
      outcome_model = outcome_model, estimand = estimand,
      control_group = control_group, doubly_robust = doubly_robust,
      nboot = nboot, boot_type = boot_type, alpha = alpha,
      anticipation = anticipation, parallel = parallel, pl_cores = pl_cores,
      data_type = data_type, cluster_var = cluster_var
    )
  } else {
    attgt_df <- .sandwich_se(
      attgt_df = attgt_df, data = data, yname = yname, tname = tname,
      idname = idname, gname = gname, xformla = xformla,
      outcome_model = outcome_model, alpha = alpha,
      data_type = data_type, cluster_var = cluster_var
    )
  }

  out <- list(
    attgt = attgt_df,
    call  = match.call(),
    args  = list(
      yname         = yname,
      tname         = tname,
      idname        = idname,
      gname         = gname,
      xformla       = xformla,
      outcome_model = outcome_model,
      estimand      = estimand,
      control_group = control_group,
      data_type     = data_type,
      weightsname   = weightsname,
      cluster_var   = cluster_var,
      doubly_robust = doubly_robust,
      alpha         = alpha,
      tlist         = tlist,
      glist         = glist
    )
  )
  class(out) <- "nonlinear_attgt"
  out
}


# ============================================================
# Dispatch: single (g,t) computation
# ============================================================

.compute_attgt_single <- function(data, yname, tname, idname, gname,
                                  xformla, g, t, base_period,
                                  outcome_model, estimand, control_group,
                                  doubly_robust, data_type) {

  if (data_type == "repeated_cross_section") {
    return(.compute_attgt_rcs(
      data = data, yname = yname, tname = tname, gname = gname,
      xformla = xformla, g = g, t = t, base_period = base_period,
      outcome_model = outcome_model, estimand = estimand,
      control_group = control_group, doubly_robust = doubly_robust
    ))
  }

  # ---- Panel path (preserved from v0.1.0; weights honoured) ----
  treated_ids <- data[[idname]][data[[gname]] == g]

  if (control_group == "nevertreated") {
    control_ids <- data[[idname]][data[[gname]] == 0]
  } else {
    control_ids <- data[[idname]][data[[gname]] == 0 | data[[gname]] > t]
    control_ids <- setdiff(control_ids, treated_ids)
  }

  if (length(treated_ids) == 0 || length(control_ids) == 0) {
    return(list(g = g, t = t, att = NA_real_, converged = FALSE,
                n_treated = 0L, n_control = 0L))
  }

  relevant_ids     <- c(treated_ids, control_ids)
  relevant_periods <- c(base_period, t)

  sub <- data[data[[idname]] %in% relevant_ids &
                data[[tname]]  %in% relevant_periods, , drop = FALSE]

  if (nrow(sub) == 0) {
    return(list(g = g, t = t, att = NA_real_, converged = FALSE,
                n_treated = length(treated_ids), n_control = length(control_ids)))
  }

  sub$D    <- as.integer(sub[[gname]] == g)
  sub$post <- as.integer(sub[[tname]] == t)

  sub_wide <- .make_wide(sub, idname = idname, tname = tname,
                         yname = yname, t0 = base_period, t1 = t)
  sub_wide$D <- as.integer(sub_wide[[gname]] == g)

  if (nrow(sub_wide) < 4) {
    return(list(g = g, t = t, att = NA_real_, converged = FALSE,
                n_treated = sum(sub_wide$D), n_control = sum(1 - sub_wide$D)))
  }

  result <- tryCatch({
    if (doubly_robust) {
      .dr_attgt(sub_wide = sub_wide, yname = yname, xformla = xformla,
                outcome_model = outcome_model, estimand = estimand)
    } else {
      .or_attgt(sub_wide = sub_wide, yname = yname, xformla = xformla,
                outcome_model = outcome_model, estimand = estimand)
    }
  }, error = function(e) {
    list(att = NA_real_, converged = FALSE, msg = conditionMessage(e))
  })

  list(
    g         = g,
    t         = t,
    att       = result$att,
    converged = isTRUE(result$converged),
    n_treated = sum(sub_wide$D == 1),
    n_control = sum(sub_wide$D == 0)
  )
}


# ============================================================
# Repeated Cross-Section ATT(g,t)
# ============================================================

#' @keywords internal
.compute_attgt_rcs <- function(data, yname, tname, gname, xformla,
                                g, t, base_period, outcome_model,
                                estimand, control_group, doubly_robust) {

  if (control_group == "nevertreated") {
    sub <- data[(data[[gname]] == g | data[[gname]] == 0) &
                  data[[tname]] %in% c(base_period, t), , drop = FALSE]
  } else {
    ctrl_flag <- data[[gname]] == 0 | (data[[gname]] > t & data[[gname]] != g)
    sub <- data[(data[[gname]] == g | ctrl_flag) &
                  data[[tname]] %in% c(base_period, t), , drop = FALSE]
  }

  if (nrow(sub) == 0) {
    return(list(g = g, t = t, att = NA_real_, converged = FALSE,
                n_treated = 0L, n_control = 0L))
  }

  sub$D    <- as.integer(sub[[gname]] == g)
  sub$post <- as.integer(sub[[tname]] == t)
  sub$did  <- sub$D * sub$post

  n_treat <- sum(sub$D == 1 & sub$post == 1)
  n_ctrl  <- sum(sub$D == 0 & sub$post == 1)

  cell_min <- sapply(list(
    sub$D == 1 & sub$post == 1,
    sub$D == 1 & sub$post == 0,
    sub$D == 0 & sub$post == 1,
    sub$D == 0 & sub$post == 0
  ), sum)

  if (any(cell_min < 2)) {
    return(list(g = g, t = t, att = NA_real_, converged = FALSE,
                n_treated = n_treat, n_control = n_ctrl))
  }

  fam   <- .outcome_family(outcome_model)
  xvars <- labels(stats::terms(xformla))
  sw    <- sub$.sw

  result <- tryCatch({
    if (!doubly_robust) {
      .rcs_pooled_qmle(sub, yname, xvars, fam, outcome_model, estimand, sw = sw)
    } else {
      .rcs_ipw_qmle(sub, yname, xvars, xformla, fam, outcome_model, estimand, sw = sw)
    }
  }, error = function(e) {
    list(att = NA_real_, converged = FALSE)
  })

  list(
    g         = g,
    t         = t,
    att       = result$att,
    converged = isTRUE(result$converged),
    n_treated = n_treat,
    n_control = n_ctrl
  )
}


# ============================================================
# RCS Estimators
# ============================================================

# Non-DR: Wooldridge (2023) pooled QMLE with D*post interaction.
.rcs_pooled_qmle <- function(sub, yname, xvars, fam, outcome_model, estimand, sw) {
  fmla <- .rcs_formula(yname, xvars)
  environment(fmla) <- environment()
  fit  <- suppressWarnings(stats::glm(fmla, data = sub, family = fam,
                                       weights = sw))
  att  <- .rcs_extract_att(fit, sub, outcome_model, estimand, fam)
  list(att = att, converged = isTRUE(fit$converged))
}


# DR: IPW-augmented pooled QMLE.
# - propensity score P(D=1|X) fit with sampling weights `sw`
# - control IPW factor e/(1-e) normalised so weighted control mass equals
#   weighted treated mass
# - final GLM uses sw * ipw_factor
.rcs_ipw_qmle <- function(sub, yname, xvars, xformla, fam, outcome_model, estimand, sw) {

  D <- sub$D

  ps_fmla <- stats::update(xformla, D ~ .)
  environment(ps_fmla) <- environment()
  ps_fit  <- suppressWarnings(
    stats::glm(ps_fmla, data = sub, family = stats::binomial(link = "logit"),
               weights = sw))
  ehat    <- pmin(pmax(stats::predict(ps_fit, type = "response"), 1e-6),
                  1 - 1e-6)

  ipw_raw <- ifelse(D == 1, 1.0, ehat / (1.0 - ehat))
  w_D1    <- sum(sw[D == 1])
  w_D0_w  <- sum((sw * ipw_raw)[D == 0])
  if (w_D0_w < 1e-10) {
    return(.rcs_pooled_qmle(sub, yname, xvars, fam, outcome_model, estimand, sw))
  }
  ipw_factor  <- ifelse(D == 1, 1.0, ipw_raw * w_D1 / w_D0_w)
  combined_wt <- sw * ipw_factor

  fmla <- .rcs_formula(yname, xvars)
  environment(fmla) <- environment()
  fit  <- suppressWarnings(
    stats::glm(fmla, data = sub, family = fam, weights = combined_wt))
  att  <- .rcs_extract_att(fit, sub, outcome_model, estimand, fam)
  list(att = att, converged = isTRUE(fit$converged))
}


.rcs_formula <- function(yname, xvars) {
  base <- "D + post + did"
  rhs  <- if (length(xvars) > 0)
    paste(base, "+", paste(xvars, collapse = " + ")) else base
  stats::as.formula(paste(yname, "~", rhs))
}


.rcs_extract_att <- function(fit, sub, outcome_model, estimand, fam) {
  coefs <- stats::coef(fit)
  if (!("did" %in% names(coefs)) || is.na(coefs["did"])) return(NA_real_)
  att_link <- unname(coefs["did"])

  if (estimand == "ape") {
    dlink <- fam$mu.eta
    lp    <- stats::predict(fit, type = "link")
    att   <- att_link * mean(dlink(lp), na.rm = TRUE)
  } else if (estimand == "odds_ratio") {
    att <- exp(att_link)
  } else {
    att <- att_link
  }
  att
}


# ============================================================
# Shared helper: GLM family from outcome_model string
# ============================================================

#' @keywords internal
.outcome_family <- function(outcome_model) {
  switch(outcome_model,
    logit   = stats::binomial(link = "logit"),
    probit  = stats::binomial(link = "probit"),
    poisson = stats::poisson(link = "log"),
    negbin  = MASS::negative.binomial(theta = 1),
    linear  = stats::gaussian(),
    stats::gaussian()
  )
}


# ============================================================
# Panel DR estimators — v0.1.0 logic, augmented with sampling weights
# ============================================================

.dr_attgt <- function(sub_wide, yname, xformla, outcome_model, estimand) {

  Y0 <- sub_wide[[paste0(yname, "_t0")]]
  Y1 <- sub_wide[[paste0(yname, "_t1")]]
  D  <- sub_wide$D
  sw <- if (".sw" %in% names(sub_wide)) sub_wide$.sw else rep(1.0, nrow(sub_wide))

  ps_formula <- stats::update(xformla, D ~ .)
  environment(ps_formula) <- environment()
  ps_fit     <- suppressWarnings(
    stats::glm(ps_formula, data = sub_wide,
               family = stats::binomial(link = "logit"),
               weights = sw))
  pscore <- stats::predict(ps_fit, type = "response")
  pscore <- pmin(pmax(pscore, 1e-6), 1 - 1e-6)

  att <- switch(estimand,
    "att" = {
      .dr_att_linear(Y0, Y1, D, pscore, sw, sub_wide, xformla, outcome_model)
    },
    "odds_ratio" = {
      .dr_att_oddsratio(Y0, Y1, D, pscore, sw, sub_wide, xformla, outcome_model)
    },
    "ape" = {
      .dr_att_ape(Y0, Y1, D, pscore, sw, sub_wide, xformla, outcome_model)
    }
  )

  list(att = att, converged = TRUE)
}


.dr_att_linear <- function(Y0, Y1, D, pscore, sw, sub_wide, xformla, outcome_model) {

  DeltaY <- Y1 - Y0

  controls_df       <- sub_wide[D == 0, , drop = FALSE]
  controls_df$DY    <- Y1[D == 0] - Y0[D == 0]
  sw_c              <- sw[D == 0]
  or_formula        <- stats::update(xformla, DY ~ .)
  environment(or_formula) <- environment()
  or_fit            <- suppressWarnings(
    stats::lm(or_formula, data = controls_df, weights = sw_c))
  mu_hat            <- stats::predict(or_fit, newdata = sub_wide)

  pD_w     <- stats::weighted.mean(D, sw)
  w_treat  <- sw * D / pD_w
  w_cont   <- sw * (1 - D) * pscore / (1 - pscore) / pD_w

  dr_score <- w_treat * DeltaY - w_cont * DeltaY + (w_treat - w_cont) * mu_hat
  sum(dr_score) / sum(sw)
}


.dr_att_oddsratio <- function(Y0, Y1, D, pscore, sw, sub_wide, xformla, outcome_model) {

  eps <- 1e-6
  Y0_t <- pmin(pmax(Y0, eps), 1 - eps)
  Y1_t <- pmin(pmax(Y1, eps), 1 - eps)

  delta_logit <- stats::qlogis(Y1_t) - stats::qlogis(Y0_t)

  pD_w    <- stats::weighted.mean(D, sw)
  w_treat <- sw * D / pD_w
  w_cont  <- sw * (1 - D) * pscore / (1 - pscore) / pD_w

  controls_df         <- sub_wide[D == 0, , drop = FALSE]
  controls_df$delta_l <- delta_logit[D == 0]
  sw_c                <- sw[D == 0]
  or_formula          <- stats::update(xformla, delta_l ~ .)
  environment(or_formula) <- environment()
  or_fit              <- suppressWarnings(
    stats::lm(or_formula, data = controls_df, weights = sw_c))
  mu_hat              <- stats::predict(or_fit, newdata = sub_wide)

  dr_score <- w_treat * delta_logit - w_cont * delta_logit +
    (w_treat - w_cont) * mu_hat
  exp(sum(dr_score) / sum(sw))
}


.dr_att_ape <- function(Y0, Y1, D, pscore, sw, sub_wide, xformla, outcome_model) {

  if (outcome_model %in% c("logit", "probit")) {
    controls_df         <- sub_wide[D == 0, , drop = FALSE]
    controls_df$Y0c     <- Y0[D == 0]
    controls_df$Y1c     <- Y1[D == 0]
    sw_c                <- sw[D == 0]

    f0_formula <- stats::update(xformla, Y0c ~ .)
    f1_formula <- stats::update(xformla, Y1c ~ .)
    environment(f0_formula) <- environment()
    environment(f1_formula) <- environment()

    fit0 <- suppressWarnings(stats::glm(f0_formula, data = controls_df,
                                        family = stats::binomial(),
                                        weights = sw_c))
    fit1 <- suppressWarnings(stats::glm(f1_formula, data = controls_df,
                                        family = stats::binomial(),
                                        weights = sw_c))

    treated_df <- sub_wide[D == 1, , drop = FALSE]
    sw_t       <- sw[D == 1]
    p0_cf <- stats::predict(fit0, newdata = treated_df, type = "response")
    p1_cf <- stats::predict(fit1, newdata = treated_df, type = "response")

    sum(sw_t * ((Y1[D == 1] - Y0[D == 1]) - (p1_cf - p0_cf))) / sum(sw_t)
  } else {
    .dr_att_linear(Y0, Y1, D, pscore, sw, sub_wide, xformla, outcome_model)
  }
}


.or_attgt <- function(sub_wide, yname, xformla, outcome_model, estimand) {

  Y0 <- sub_wide[[paste0(yname, "_t0")]]
  Y1 <- sub_wide[[paste0(yname, "_t1")]]
  D  <- sub_wide$D
  sw <- if (".sw" %in% names(sub_wide)) sub_wide$.sw else rep(1.0, nrow(sub_wide))

  controls_df         <- sub_wide[D == 0, , drop = FALSE]
  controls_df$DY      <- Y1[D == 0] - Y0[D == 0]
  sw_c                <- sw[D == 0]
  or_formula          <- stats::update(xformla, DY ~ .)
  environment(or_formula) <- environment()
  or_fit              <- suppressWarnings(
    stats::lm(or_formula, data = controls_df, weights = sw_c))

  treated_df  <- sub_wide[D == 1, , drop = FALSE]
  mu_hat_t    <- stats::predict(or_fit, newdata = treated_df)
  sw_t        <- sw[D == 1]

  att <- sum(sw_t * (Y1[D == 1] - Y0[D == 1] - mu_hat_t)) / sum(sw_t)
  list(att = att, converged = TRUE)
}


# ============================================================
# Helper: wide format (panel only) — carries sampling weight & cluster
# ============================================================

.make_wide <- function(sub, idname, tname, yname, t0, t1) {

  base_cols_0 <- c(idname, yname, "D")
  if (".sw" %in% names(sub)) base_cols_0 <- c(base_cols_0, ".sw")
  base_cols_1 <- c(idname, yname, "D")

  sub0 <- sub[sub[[tname]] == t0, base_cols_0, drop = FALSE]
  sub1 <- sub[sub[[tname]] == t1, base_cols_1, drop = FALSE]

  names(sub0)[names(sub0) == yname] <- paste0(yname, "_t0")
  names(sub1)[names(sub1) == yname] <- paste0(yname, "_t1")
  names(sub0)[names(sub0) == "D"]   <- "D_t0"
  names(sub1)[names(sub1) == "D"]   <- "D_t1"

  merged <- merge(sub0, sub1, by = idname)
  merged$D <- merged$D_t1

  # Other carried columns (e.g. covariates, cluster_var) from t0
  extra_cols <- setdiff(names(sub),
                        c(idname, tname, yname, "D", "post", ".sw"))
  if (length(extra_cols) > 0) {
    sub0_extra <- sub[sub[[tname]] == t0, c(idname, extra_cols), drop = FALSE]
    merged <- merge(merged, sub0_extra, by = idname, all.x = TRUE)
  }

  merged
}
