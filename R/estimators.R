# Helper: build a formula safely from a base formula + covariates
.make_formula <- function(lhs_str, base_rhs_str, xformla) {
  xvars <- labels(stats::terms(xformla))
  if (length(xvars) > 0 && !(length(xvars) == 1 && xvars == "")) {
    rhs <- paste(c(base_rhs_str, xvars), collapse = " + ")
  } else {
    rhs <- base_rhs_str
  }
  stats::as.formula(paste(lhs_str, "~", rhs))
}

#' @title Binary Outcome DiD: Logit Estimator
#'
#' @description
#' Estimates a 2x2 difference-in-differences model with a binary outcome using
#' logistic regression on the log-odds scale, reporting both the log-odds DiD
#' coefficient and the average partial effect (APE) on the probability scale.
#'
#' @param data A data frame (long format).
#' @param yname Character. Binary outcome variable name.
#' @param tname Character. Time period variable name.
#' @param idname Character. Unit ID variable name.
#' @param treat_period Numeric. The treatment (post) period.
#' @param control_period Numeric. The pre-treatment baseline period.
#' @param dname Character. Treatment indicator variable name (optional).
#' @param gname Character. Cohort variable name (optional).
#' @param xformla One-sided formula for covariates. Default \code{~1}.
#' @param se_type Character. SE type: \code{"robust"} (default),
#'   \code{"cluster"}, or \code{"analytical"}.
#' @param cluster_var Character. Clustering variable (if \code{se_type = "cluster"}).
#'
#' @return A list of class \code{binary_did_logit}.
#' @examples
#' dat <- sim_binary_panel(n = 500, nperiods = 4, prop_treated = 0.5)
#' dat2 <- dat[dat$period %in% c(2, 3), ]
#' res <- binary_did_logit(dat2, yname = "y", tname = "period",
#'                          idname = "id", treat_period = 3,
#'                          control_period = 2, gname = "g")
#' print(res)
#' @export
binary_did_logit <- function(
    data, yname, tname, idname, treat_period, control_period,
    dname = NULL, gname = NULL, xformla = ~1,
    se_type = c("robust", "cluster", "analytical"), cluster_var = NULL
) {
  se_type <- match.arg(se_type)
  .binary_did_impl(data, yname, tname, idname, treat_period, control_period,
                    dname, gname, xformla, se_type, cluster_var, link = "logit")
}

#' @title Binary Outcome DiD: Probit Estimator
#'
#' @description
#' Estimates 2x2 DiD with binary outcome using probit regression.
#' Parallel trends assumed on the probit (inverse-normal) scale.
#'
#' @inheritParams binary_did_logit
#' @return A list of class \code{binary_did_probit}.
#' @examples
#' dat <- sim_binary_panel(n = 500, nperiods = 4, prop_treated = 0.5)
#' dat2 <- dat[dat$period %in% c(2, 3), ]
#' res <- binary_did_probit(dat2, "y", "period", "id", 3, 2, gname = "g")
#' print(res)
#' @export
binary_did_probit <- function(
    data, yname, tname, idname, treat_period, control_period,
    dname = NULL, gname = NULL, xformla = ~1,
    se_type = c("robust", "cluster", "analytical"), cluster_var = NULL
) {
  se_type <- match.arg(se_type)
  .binary_did_impl(data, yname, tname, idname, treat_period, control_period,
                    dname, gname, xformla, se_type, cluster_var, link = "probit")
}

#' @title Doubly-Robust Binary DiD
#'
#' @description
#' Doubly-robust estimator for binary outcomes combining a nonlinear outcome
#' regression model with inverse probability weighting via propensity score.
#' Consistent if EITHER the outcome model OR the propensity score is correctly
#' specified.
#'
#' @inheritParams binary_did_logit
#' @param outcome_model Character. \code{"logit"} (default) or \code{"probit"}.
#' @return A list of class \code{binary_did_dr}.
#' @examples
#' dat <- sim_binary_panel(n = 500, nperiods = 4, prop_treated = 0.5)
#' dat2 <- dat[dat$period %in% c(2, 3), ]
#' res <- binary_did_dr(dat2, "y", "period", "id", 3, 2, gname = "g",
#'                       outcome_model = "logit")
#' print(res)
#' @export
binary_did_dr <- function(
    data, yname, tname, idname, treat_period, control_period,
    dname = NULL, gname = NULL, xformla = ~1,
    outcome_model = c("logit", "probit"),
    se_type = c("robust", "cluster", "analytical"), cluster_var = NULL
) {
  outcome_model <- match.arg(outcome_model)
  se_type       <- match.arg(se_type)

  sub <- data[data[[tname]] %in% c(control_period, treat_period), , drop = FALSE]
  sub$D <- if (!is.null(gname)) as.integer(sub[[gname]] != 0) else
            if (!is.null(dname)) as.integer(sub[[dname]]) else
            stop("Must provide either 'gname' or 'dname'.")

  sub_wide <- .make_wide(sub, idname, tname, yname, control_period, treat_period)
  sub_wide$D <- as.integer(sub_wide$D_t1)

  Y0 <- sub_wide[[paste0(yname, "_t0")]]
  Y1 <- sub_wide[[paste0(yname, "_t1")]]
  D  <- sub_wide$D

  ps_formula <- stats::update(xformla, D ~ .)
  ps_fit     <- suppressWarnings(
    stats::glm(ps_formula, data = sub_wide, family = stats::binomial()))
  pscore <- pmin(pmax(stats::predict(ps_fit, type = "response"), 1e-6), 1 - 1e-6)

  att_val <- .dr_att_linear(Y0, Y1, D, pscore, sub_wide, xformla, outcome_model)

  DeltaY <- Y1 - Y0
  p_D    <- mean(D)
  w_t    <- D / p_D
  w_c    <- (1 - D) * pscore / (1 - pscore) / p_D

  ctrl_df     <- sub_wide[D == 0, , drop = FALSE]
  ctrl_df$DY  <- DeltaY[D == 0]
  or_fit      <- suppressWarnings(stats::lm(.make_formula("DY", "1", xformla),
                                             data = ctrl_df))
  mu_hat <- stats::predict(or_fit, newdata = sub_wide)
  psi    <- w_t * (DeltaY - mu_hat) - w_c * (DeltaY - mu_hat)
  se_val <- stats::sd(psi) / sqrt(length(psi))

  out <- list(
    att = att_val, se = se_val,
    ci_lo = att_val - 1.96 * se_val, ci_hi = att_val + 1.96 * se_val,
    tstat = att_val / se_val,
    pval  = 2 * stats::pnorm(-abs(att_val / se_val)),
    n = nrow(sub_wide), n_treat = sum(D), n_cont = sum(1 - D),
    outcome_model = outcome_model
  )
  class(out) <- "binary_did_dr"
  out
}

#' @title Count Outcome DiD: Poisson Estimator
#'
#' @description
#' Estimates DiD for count outcomes using a Poisson quasi-maximum likelihood
#' (QMLE) estimator with a log-linear parallel trends assumption.
#' The treatment effect is a multiplicative rate ratio.
#'
#' @inheritParams binary_did_logit
#' @param offset Character. Name of offset variable. Default \code{NULL}.
#' @return A list of class \code{count_did_poisson}.
#' @examples
#' dat <- sim_count_panel(n = 400, nperiods = 6, prop_treated = 0.4)
#' dat2 <- dat[dat$period %in% c(2, 4), ]
#' res <- count_did_poisson(dat2, "y", "period", "id", 4, 2, gname = "g")
#' print(res)
#' @export
count_did_poisson <- function(
    data, yname, tname, idname, treat_period, control_period,
    dname = NULL, gname = NULL, xformla = ~1, offset = NULL,
    se_type = c("robust", "cluster", "analytical"), cluster_var = NULL
) {
  se_type <- match.arg(se_type)

  sub <- data[data[[tname]] %in% c(control_period, treat_period), , drop = FALSE]
  sub$D <- if (!is.null(gname)) as.integer(sub[[gname]] != 0) else
            if (!is.null(dname)) as.integer(sub[[dname]]) else
            stop("Must provide 'gname' or 'dname'.")

  sub$post <- as.integer(sub[[tname]] == treat_period)
  sub$did  <- sub$D * sub$post

  fit_formula <- .make_formula(yname, "D + post + did", xformla)

  if (!is.null(offset)) {
    fit <- suppressWarnings(
      stats::glm(fit_formula, data = sub, family = stats::poisson(link = "log"),
                 offset = log(sub[[offset]] + 1))
    )
  } else {
    fit <- suppressWarnings(
      stats::glm(fit_formula, data = sub, family = stats::quasipoisson(link = "log"))
    )
  }

  vcov_mat <- if (se_type == "robust") sandwich::vcovHC(fit, type = "HC1") else stats::vcov(fit)
  coef_tab <- lmtest::coeftest(fit, vcov = vcov_mat)
  did_coef <- stats::coef(fit)["did"]
  did_se   <- sqrt(vcov_mat["did", "did"])
  rate_ratio <- exp(did_coef)
  mean_y0    <- mean(sub[[yname]][sub$D == 1 & sub$post == 0], na.rm = TRUE)

  out <- list(
    att_log_rr = did_coef, se_log_rr = did_se,
    rate_ratio = rate_ratio,
    att_ape    = mean_y0 * (rate_ratio - 1),
    ci_lo_log  = did_coef - 1.96 * did_se,
    ci_hi_log  = did_coef + 1.96 * did_se,
    ci_lo_rr   = exp(did_coef - 1.96 * did_se),
    ci_hi_rr   = exp(did_coef + 1.96 * did_se),
    tstat      = did_coef / did_se,
    pval       = 2 * stats::pnorm(-abs(did_coef / did_se)),
    fit        = fit, coef_table = coef_tab
  )
  class(out) <- "count_did_poisson"
  out
}

#' @title Odds-Ratio DiD for Binary Outcomes
#'
#' @description
#' Estimates the odds-ratio difference-in-differences (OR-DiD) for binary
#' outcomes. OR-DiD equals 1 under no treatment effect and is invariant to
#' which group is labelled treatment.
#'
#' @inheritParams binary_did_logit
#' @return A list of class \code{odds_ratio_did}.
#' @examples
#' dat <- sim_binary_panel(n = 500, nperiods = 4, prop_treated = 0.5)
#' dat2 <- dat[dat$period %in% c(2, 3), ]
#' res <- odds_ratio_did(dat2, "y", "period", "id", 3, 2, gname = "g")
#' print(res)
#' @export
odds_ratio_did <- function(
    data, yname, tname, idname, treat_period, control_period,
    dname = NULL, gname = NULL, xformla = ~1
) {
  sub <- data[data[[tname]] %in% c(control_period, treat_period), , drop = FALSE]
  sub$D <- if (!is.null(gname)) as.integer(sub[[gname]] != 0) else
            if (!is.null(dname)) as.integer(sub[[dname]]) else
            stop("Must provide 'gname' or 'dname'.")

  sub$post <- as.integer(sub[[tname]] == treat_period)
  sub$did  <- sub$D * sub$post

  fit_formula <- .make_formula(yname, "D + post + did", xformla)
  fit <- suppressWarnings(
    stats::glm(fit_formula, data = sub, family = stats::binomial(link = "logit"))
  )

  vcov_rob   <- sandwich::vcovHC(fit, type = "HC1")
  log_or_did <- stats::coef(fit)["did"]
  se_log_or  <- sqrt(vcov_rob["did", "did"])
  or_did     <- exp(log_or_did)

  out <- list(
    log_or_did = log_or_did, se_log_or = se_log_or,
    or_did     = or_did,
    ci_lo_log  = log_or_did - 1.96 * se_log_or,
    ci_hi_log  = log_or_did + 1.96 * se_log_or,
    ci_lo_or   = exp(log_or_did - 1.96 * se_log_or),
    ci_hi_or   = exp(log_or_did + 1.96 * se_log_or),
    tstat      = log_or_did / se_log_or,
    pval       = 2 * stats::pnorm(-abs(log_or_did / se_log_or)),
    fit        = fit,
    interpretation = paste0(
      "OR-DiD = ", round(or_did, 3),
      " [95% CI: ", round(exp(log_or_did - 1.96 * se_log_or), 3),
      ", ", round(exp(log_or_did + 1.96 * se_log_or), 3), "]"
    )
  )
  class(out) <- "odds_ratio_did"
  out
}

# ============================================================
# Internal implementation
# ============================================================

.binary_did_impl <- function(data, yname, tname, idname, treat_period,
                               control_period, dname, gname, xformla,
                               se_type, cluster_var, link) {

  sub <- data[data[[tname]] %in% c(control_period, treat_period), , drop = FALSE]
  sub$D <- if (!is.null(gname)) as.integer(sub[[gname]] != 0) else
            if (!is.null(dname)) as.integer(sub[[dname]]) else
            stop("Must provide either 'gname' or 'dname'.")

  sub$post <- as.integer(sub[[tname]] == treat_period)
  sub$did  <- sub$D * sub$post

  fit_formula <- .make_formula(yname, "D + post + did", xformla)
  fit <- suppressWarnings(
    stats::glm(fit_formula, data = sub, family = stats::binomial(link = link))
  )

  vcov_mat <- if (se_type == "robust") {
    sandwich::vcovHC(fit, type = "HC1")
  } else if (se_type == "cluster" && !is.null(cluster_var)) {
    sandwich::vcovCL(fit, cluster = sub[[cluster_var]])
  } else {
    stats::vcov(fit)
  }

  coef_tab <- lmtest::coeftest(fit, vcov = vcov_mat)
  did_coef <- stats::coef(fit)["did"]
  did_se   <- sqrt(vcov_mat["did", "did"])

  link_fn  <- if (link == "logit") stats::plogis else stats::pnorm
  dlink_fn <- if (link == "logit") {
    function(x) stats::plogis(x) * (1 - stats::plogis(x))
  } else {
    stats::dnorm
  }
  lp       <- stats::predict(fit, type = "link")
  ape_mult <- mean(dlink_fn(lp))
  att_ape  <- did_coef * ape_mult

  out <- list(
    att_link   = did_coef, se_link = did_se,
    att_ape    = att_ape,
    ci_lo_link = did_coef - 1.96 * did_se,
    ci_hi_link = did_coef + 1.96 * did_se,
    tstat = did_coef / did_se,
    pval  = 2 * stats::pnorm(-abs(did_coef / did_se)),
    fit = fit, coef_table = coef_tab, link = link
  )
  class(out) <- paste0("binary_did_", link)
  out
}
