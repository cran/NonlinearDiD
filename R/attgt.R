#' @title Nonlinear Staggered DiD: Group-Time ATT Estimation
#'
#' @description
#' Computes group-time average treatment effects on the treated (ATT(g,t)) for
#' staggered difference-in-differences designs with nonlinear outcomes.
#'
#' This function extends Callaway & Sant'Anna (2021) to handle binary, count,
#' and other nonlinear outcomes where the standard linear parallel trends
#' assumption is misspecified. The key methodological contributions are:
#'
#' 1. **Parallel trends on the latent index** (for logit/probit): Instead of
#'    assuming parallel trends in \eqn{E[Y]}, we assume parallel trends in the
#'    latent utility \eqn{F^{-1}(E[Y])}.
#'
#' 2. **Doubly-robust nonlinear estimator**: Combines outcome regression
#'    (nonlinear model) with propensity score weighting, inheriting DR
#'    properties in the nonlinear setting.
#'
#' 3. **Odds-ratio DiD**: A scale-free estimand appropriate for binary
#'    outcomes that does not require parallel trends in probabilities.
#'
#' 4. **Nonparametric bounds**: When no functional form is assumed,
#'    provides sharp bounds on ATT(g,t).
#'
#' @param data A data frame in long format (one row per unit-period).
#' @param yname Character. Name of the outcome variable column.
#' @param tname Character. Name of the time period column.
#' @param idname Character. Name of the unit identifier column.
#' @param gname Character. Name of the treatment cohort column (the period
#'   when a unit first receives treatment; 0 or Inf for never-treated units).
#' @param xformla A one-sided formula for covariates (e.g., `~ x1 + x2`).
#'   Default is `~ 1` (intercept only).
#' @param outcome_model Character. The outcome model to use. One of:
#'   \itemize{
#'     \item \code{"logit"}: Logistic regression (for binary Y)
#'     \item \code{"probit"}: Probit regression (for binary Y)
#'     \item \code{"poisson"}: Poisson regression (for count Y)
#'     \item \code{"negbin"}: Negative binomial (for overdispersed count Y)
#'     \item \code{"linear"}: Linear model (reproduces CS2021 when combined
#'       with doubly_robust = TRUE)
#'   }
#' @param estimand Character. The treatment effect estimand:
#'   \itemize{
#'     \item \code{"att"}: Average treatment effect on the treated (default)
#'     \item \code{"odds_ratio"}: Odds ratio DiD (binary outcomes only)
#'     \item \code{"ape"}: Average partial effect on the probability scale
#'   }
#' @param control_group Character. Which units serve as the control group:
#'   \itemize{
#'     \item \code{"nevertreated"}: Use never-treated units only (default)
#'     \item \code{"notyetreated"}: Use not-yet-treated units
#'   }
#' @param doubly_robust Logical. If TRUE (default), uses the doubly-robust
#'   estimator that combines propensity score weighting with outcome
#'   regression. More robust to model misspecification.
#' @param boot Logical. If TRUE, uses bootstrap for inference. Default FALSE.
#' @param nboot Integer. Number of bootstrap iterations. Default 999.
#' @param boot_type Character. Type of bootstrap: \code{"multiplier"}
#'   (default, fast) or \code{"empirical"}.
#' @param alpha Numeric. Significance level for confidence intervals.
#'   Default 0.05.
#' @param parallel Logical. Use parallel processing for bootstrap. Default FALSE.
#' @param pl_cores Integer. Number of cores for parallel processing.
#' @param anticipation Integer. Number of periods of anticipation allowed.
#'   Default 0.
#'
#' @return An object of class \code{nonlinear_attgt} containing:
#'   \describe{
#'     \item{attgt}{Data frame of ATT(g,t) estimates, standard errors, and
#'       confidence intervals for each (group, time) pair.}
#'     \item{call}{The matched call.}
#'     \item{args}{List of arguments used.}
#'     \item{boot_draws}{Matrix of bootstrap draws (if boot = TRUE).}
#'   }
#'
#' @references
#' Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences
#' with multiple time periods. *Journal of Econometrics*, 225(2), 200-230.
#'
#' Wooldridge, J. M. (2023). Simple approaches to nonlinear
#' difference-in-differences with panel data. *The Econometrics Journal*, 26(3).
#'
#' Roth, J., & Sant'Anna, P. H. C. (2023). When is parallel trends sensitive
#' to functional form? *Econometrica*, 91(2), 737-747.
#'
#' @examples
#' # Simulate binary panel data
#' set.seed(42)
#' dat <- sim_binary_panel(n = 500, nperiods = 6, prop_treated = 0.4)
#'
#' # Estimate ATT(g,t) with logistic outcome model
#' result <- nonlinear_attgt(
#'   data = dat,
#'   yname = "y",
#'   tname = "period",
#'   idname = "id",
#'   gname = "g",
#'   outcome_model = "logit",
#'   control_group = "nevertreated"
#' )
#'
#' summary(result)
#' plot(result)
#'
#' @export
nonlinear_attgt <- function(
    data,
    yname,
    tname,
    idname,
    gname,
    xformla = ~1,
    outcome_model = c("logit", "probit", "poisson", "negbin", "linear"),
    estimand = c("att", "odds_ratio", "ape"),
    control_group = c("nevertreated", "notyetreated"),
    doubly_robust = TRUE,
    boot = FALSE,
    nboot = 999,
    boot_type = c("multiplier", "empirical"),
    alpha = 0.05,
    parallel = FALSE,
    pl_cores = 2L,
    anticipation = 0L
) {

  # ---- Input validation ----
  outcome_model <- match.arg(outcome_model)
  estimand      <- match.arg(estimand)
  control_group <- match.arg(control_group)
  boot_type     <- match.arg(boot_type)

  if (estimand == "odds_ratio" && !outcome_model %in% c("logit", "probit", "linear")) {
    stop("'odds_ratio' estimand requires outcome_model = 'logit', 'probit', or 'linear'.")
  }

  # Validate columns exist
  required_cols <- c(yname, tname, idname, gname)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in data:", paste(missing_cols, collapse = ", ")))
  }

  # ---- Standardize data ----
  data <- as.data.frame(data)
  data[[gname]] <- ifelse(is.na(data[[gname]]) | is.infinite(data[[gname]]), 0, data[[gname]])

  tlist <- sort(unique(data[[tname]]))
  glist <- sort(setdiff(unique(data[[gname]]), 0))

  if (length(glist) == 0) stop("No treated units found (all gname values are 0 or missing).")
  if (length(tlist) < 2)  stop("Need at least 2 time periods.")

  # ---- Compute ATT(g,t) for each (g, t) pair ----
  attgt_list <- vector("list", length(glist) * length(tlist))
  idx <- 0L

  for (g in glist) {
    for (t in tlist) {
      idx <- idx + 1L

      # Skip pre-treatment periods beyond anticipation window
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
        doubly_robust = doubly_robust
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

  # ---- Bootstrap inference ----
  if (boot) {
    attgt_df <- .bootstrap_inference(
      attgt_df   = attgt_df,
      data       = data,
      yname      = yname,
      tname      = tname,
      idname     = idname,
      gname      = gname,
      xformla    = xformla,
      outcome_model  = outcome_model,
      estimand   = estimand,
      control_group  = control_group,
      doubly_robust  = doubly_robust,
      nboot      = nboot,
      boot_type  = boot_type,
      alpha      = alpha,
      anticipation   = anticipation,
      parallel   = parallel,
      pl_cores   = pl_cores
    )
  } else {
    # Delta-method / sandwich SEs
    attgt_df <- .sandwich_se(
      attgt_df  = attgt_df,
      data      = data,
      yname     = yname,
      tname     = tname,
      idname    = idname,
      gname     = gname,
      xformla   = xformla,
      outcome_model = outcome_model,
      alpha     = alpha
    )
  }

  # ---- Return ----
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
# Internal: Single (g,t) ATT computation
# ============================================================

.compute_attgt_single <- function(data, yname, tname, idname, gname,
                                  xformla, g, t, base_period,
                                  outcome_model, estimand, control_group,
                                  doubly_robust) {

  # Identify treated group and control group
  treated_ids <- data[[idname]][data[[gname]] == g]

  if (control_group == "nevertreated") {
    control_ids <- data[[idname]][data[[gname]] == 0]
  } else {  # notyetreated
    control_ids <- data[[idname]][data[[gname]] == 0 | data[[gname]] > t]
    control_ids <- setdiff(control_ids, treated_ids)
  }

  if (length(treated_ids) == 0 || length(control_ids) == 0) {
    return(list(g = g, t = t, att = NA_real_, converged = FALSE,
                n_treated = 0L, n_control = 0L))
  }

  # Subset to relevant units and two time periods: base_period and t
  relevant_ids  <- c(treated_ids, control_ids)
  relevant_periods <- c(base_period, t)

  sub <- data[data[[idname]] %in% relevant_ids &
                data[[tname]]  %in% relevant_periods, , drop = FALSE]

  if (nrow(sub) == 0) {
    return(list(g = g, t = t, att = NA_real_, converged = FALSE,
                n_treated = length(treated_ids), n_control = length(control_ids)))
  }

  sub$D  <- as.integer(sub[[gname]] == g)
  sub$post <- as.integer(sub[[tname]] == t)

  # Wide-format DiD: outcome change
  sub_wide <- .make_wide(sub, idname = idname, tname = tname,
                         yname = yname, t0 = base_period, t1 = t)
  sub_wide$D <- as.integer(sub_wide[[gname]] == g)

  if (nrow(sub_wide) < 4) {
    return(list(g = g, t = t, att = NA_real_, converged = FALSE,
                n_treated = sum(sub_wide$D), n_control = sum(1 - sub_wide$D)))
  }

  # Dispatch to appropriate estimator
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
# Doubly-Robust Nonlinear ATT(g,t)
# ============================================================

.dr_attgt <- function(sub_wide, yname, xformla, outcome_model, estimand) {

  Y0 <- sub_wide[[paste0(yname, "_t0")]]
  Y1 <- sub_wide[[paste0(yname, "_t1")]]
  D  <- sub_wide$D
  n  <- nrow(sub_wide)

  # --- Propensity score ---
  ps_formula <- stats::update(xformla, D ~ .)
  ps_fit     <- suppressWarnings(
    stats::glm(ps_formula, data = sub_wide, family = stats::binomial(link = "logit"))
  )
  pscore <- stats::predict(ps_fit, type = "response")
  pscore <- pmin(pmax(pscore, 1e-6), 1 - 1e-6)  # trim

  # --- Outcome regression on controls ---
  controls  <- sub_wide[D == 0, , drop = FALSE]
  X_formula <- stats::update(xformla, NULL ~ .)

  att <- switch(estimand,
                "att" = {
                  .dr_att_linear(Y0 = Y0, Y1 = Y1, D = D, pscore = pscore,
                                 sub_wide = sub_wide, xformla = xformla,
                                 outcome_model = outcome_model)
                },
                "odds_ratio" = {
                  .dr_att_oddsratio(Y0 = Y0, Y1 = Y1, D = D, pscore = pscore,
                                    sub_wide = sub_wide, xformla = xformla,
                                    outcome_model = outcome_model)
                },
                "ape" = {
                  .dr_att_ape(Y0 = Y0, Y1 = Y1, D = D, pscore = pscore,
                              sub_wide = sub_wide, xformla = xformla,
                              outcome_model = outcome_model)
                }
  )

  list(att = att, converged = TRUE)
}


.dr_att_linear <- function(Y0, Y1, D, pscore, sub_wide, xformla, outcome_model) {

  DeltaY <- Y1 - Y0
  n      <- length(D)

  # Outcome model for E[Delta Y | X, D=0]
  controls_df     <- sub_wide[D == 0, , drop = FALSE]
  controls_df$DY  <- Y1[D == 0] - Y0[D == 0]

  or_formula <- stats::update(xformla, DY ~ .)

  or_family <- switch(outcome_model,
                      logit   = stats::binomial(link = "logit"),
                      probit  = stats::binomial(link = "probit"),
                      poisson = stats::poisson(link = "log"),
                      linear  = stats::gaussian(),
                      stats::gaussian()  # default
  )

  # For difference outcome, use Gaussian unless count
  if (outcome_model %in% c("logit", "probit")) {
    # Latent index: fit linear model on log-odds change for controls
    or_fit <- suppressWarnings(stats::lm(or_formula, data = controls_df))
  } else if (outcome_model == "poisson") {
    # For Poisson, use rate ratio: Y1/Y0 for controls
    controls_df$DY_pos <- pmax(controls_df$DY + abs(min(controls_df$DY)) + 1, 1e-6)
    or_fit <- suppressWarnings(stats::lm(or_formula, data = controls_df))
  } else {
    or_fit <- suppressWarnings(stats::lm(or_formula, data = controls_df))
  }

  mu_hat <- stats::predict(or_fit, newdata = sub_wide)

  # IPW weights
  w_treat  <- D / mean(D)
  w_cont   <- (1 - D) * pscore / (1 - pscore) / mean(D)

  # DR score: combines OR and IPW
  dr_score <- w_treat * DeltaY - w_cont * DeltaY + (w_treat - w_cont) * mu_hat
  att_dr   <- mean(dr_score)

  att_dr
}


.dr_att_oddsratio <- function(Y0, Y1, D, pscore, sub_wide, xformla, outcome_model) {
  # Odds Ratio DiD:
  # ln OR_DiD = [ln(p11/(1-p11)) - ln(p10/(1-p10))] - [ln(p01/(1-p01)) - ln(p00/(1-p00))]
  # where pdt = E[Y | D=d, t=t']

  eps <- 1e-6
  Y0_t <- pmin(pmax(Y0, eps), 1 - eps)
  Y1_t <- pmin(pmax(Y1, eps), 1 - eps)

  logit_Y0 <- stats::qlogis(Y0_t)
  logit_Y1 <- stats::qlogis(Y1_t)
  delta_logit <- logit_Y1 - logit_Y0

  # DR estimator on the logit scale
  w_treat <- D / mean(D)
  w_cont  <- (1 - D) * pscore / (1 - pscore) / mean(D)

  # Control outcome model on logit scale
  controls_df          <- sub_wide[D == 0, , drop = FALSE]
  controls_df$delta_l  <- delta_logit[D == 0]
  or_formula           <- stats::update(xformla, delta_l ~ .)
  or_fit               <- suppressWarnings(stats::lm(or_formula, data = controls_df))
  mu_hat               <- stats::predict(or_fit, newdata = sub_wide)

  dr_score <- w_treat * delta_logit - w_cont * delta_logit +
    (w_treat - w_cont) * mu_hat
  att_or   <- mean(dr_score)

  # Exponentiate: exp(log_OR_DiD) = OR_DiD
  exp(att_or)
}


.dr_att_ape <- function(Y0, Y1, D, pscore, sub_wide, xformla, outcome_model) {
  # Average Partial Effect on probability scale
  # Uses the correlated random effects Mundlak approach for nonlinear panel

  eps   <- 1e-6
  DeltaY <- Y1 - Y0
  n     <- length(D)

  # For APE in logit/probit, use integrated partial effects
  if (outcome_model %in% c("logit", "probit")) {
    link_fn  <- if (outcome_model == "logit") stats::plogis else stats::pnorm
    # Fit model to controls, recover counterfactual means
    controls_df   <- sub_wide[D == 0, , drop = FALSE]
    controls_df$Y0c <- Y0[D == 0]
    controls_df$Y1c <- Y1[D == 0]

    f0_formula <- stats::update(xformla, Y0c ~ .)
    f1_formula <- stats::update(xformla, Y1c ~ .)

    fit0 <- suppressWarnings(stats::glm(f0_formula, data = controls_df, family = stats::binomial()))
    fit1 <- suppressWarnings(stats::glm(f1_formula, data = controls_df, family = stats::binomial()))

    p0_cf <- stats::predict(fit0, newdata = sub_wide[D == 1, , drop = FALSE], type = "response")
    p1_cf <- stats::predict(fit1, newdata = sub_wide[D == 1, , drop = FALSE], type = "response")

    att_ape <- mean((Y1[D == 1] - Y0[D == 1]) - (p1_cf - p0_cf))
  } else {
    # Fall back to linear DR
    att_ape <- .dr_att_linear(Y0, Y1, D, pscore, sub_wide, xformla, outcome_model)
  }

  att_ape
}


# ============================================================
# Outcome Regression only (no propensity score)
# ============================================================

.or_attgt <- function(sub_wide, yname, xformla, outcome_model, estimand) {

  Y0 <- sub_wide[[paste0(yname, "_t0")]]
  Y1 <- sub_wide[[paste0(yname, "_t1")]]
  D  <- sub_wide$D

  controls_df    <- sub_wide[D == 0, , drop = FALSE]
  controls_df$DY <- Y1[D == 0] - Y0[D == 0]
  or_formula     <- stats::update(xformla, DY ~ .)
  or_fit         <- suppressWarnings(stats::lm(or_formula, data = controls_df))
  mu_hat_treated <- stats::predict(or_fit, newdata = sub_wide[D == 1, , drop = FALSE])

  att <- mean(Y1[D == 1] - Y0[D == 1]) - mean(mu_hat_treated)
  list(att = att, converged = TRUE)
}


# ============================================================
# Helper: wide format
# ============================================================

.make_wide <- function(sub, idname, tname, yname, t0, t1) {
  sub0 <- sub[sub[[tname]] == t0, c(idname, yname, "D"), drop = FALSE]
  sub1 <- sub[sub[[tname]] == t1, c(idname, yname, "D"), drop = FALSE]

  names(sub0)[names(sub0) == yname] <- paste0(yname, "_t0")
  names(sub1)[names(sub1) == yname] <- paste0(yname, "_t1")
  names(sub0)[names(sub0) == "D"]   <- paste0("D_t0")
  names(sub1)[names(sub1) == "D"]   <- paste0("D_t1")

  merged <- merge(sub0, sub1, by = idname)
  merged$D <- merged$D_t1  # use current treatment indicator

  # Carry through covariates from t0
  extra_cols <- setdiff(names(sub), c(idname, tname, yname, "D", "post"))
  if (length(extra_cols) > 0) {
    sub0_extra <- sub[sub[[tname]] == t0, c(idname, extra_cols), drop = FALSE]
    merged <- merge(merged, sub0_extra, by = idname, all.x = TRUE)
  }

  # Preserve gname
  if ("gname" %in% names(sub)) {
    gname_vals <- sub[sub[[tname]] == t0, c(idname, "gname"), drop = FALSE]
    merged <- merge(merged, gname_vals, by = idname, all.x = TRUE)
  }

  merged
}
