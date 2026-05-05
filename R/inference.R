#' @title Inference for Nonlinear DiD
#' @description Internal functions for bootstrap and delta-method
#' standard errors in nonlinear staggered DiD.
#' @keywords internal

# ============================================================
# Multiplier / Empirical Bootstrap
# ============================================================

.bootstrap_inference <- function(attgt_df, data, yname, tname, idname, gname,
                                 xformla, outcome_model, estimand, control_group,
                                 doubly_robust, nboot, boot_type, alpha,
                                 anticipation, parallel, pl_cores) {

  ids      <- unique(data[[idname]])
  n_units  <- length(ids)
  tlist    <- sort(unique(data[[tname]]))
  glist    <- sort(setdiff(unique(data[[gname]]), 0))

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(pl_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, varlist = ls(envir = environment()),
                            envir = environment())
    boot_fn <- function(b) parallel::clusterCall(cl, .single_boot,
                                                 b, ids, n_units, data,
                                                 yname, tname, idname, gname,
                                                 xformla, outcome_model,
                                                 estimand, control_group,
                                                 doubly_robust, boot_type,
                                                 anticipation, glist, tlist)
    boot_draws <- parallel::parSapply(cl, seq_len(nboot), function(b)
      .single_boot(b, ids, n_units, data, yname, tname, idname, gname,
                   xformla, outcome_model, estimand, control_group,
                   doubly_robust, boot_type, anticipation, glist, tlist))
  } else {
    boot_draws <- matrix(NA_real_, nrow = nrow(attgt_df), ncol = nboot)
    for (b in seq_len(nboot)) {
      bv <- .single_boot(b, ids, n_units, data, yname, tname, idname, gname,
                         xformla, outcome_model, estimand, control_group,
                         doubly_robust, boot_type, anticipation, glist, tlist)
      boot_draws[, b] <- bv
    }
  }

  # SE and CI from bootstrap distribution
  attgt_df$se    <- apply(boot_draws, 1, function(x) stats::sd(x, na.rm = TRUE))
  attgt_df$ci_lo <- attgt_df$att - stats::qnorm(1 - alpha / 2) * attgt_df$se
  attgt_df$ci_hi <- attgt_df$att + stats::qnorm(1 - alpha / 2) * attgt_df$se

  # Simultaneous confidence bands via Holm-Bonferroni
  attgt_df$ci_lo_simult <- attgt_df$att -
    apply(boot_draws, 1, function(x) stats::quantile(abs(x - attgt_df$att[1]),
                                                     1 - alpha, na.rm = TRUE))
  attgt_df$ci_hi_simult <- attgt_df$att +
    apply(boot_draws, 1, function(x) stats::quantile(abs(x - attgt_df$att[1]),
                                                     1 - alpha, na.rm = TRUE))

  attgt_df
}


.single_boot <- function(b, ids, n_units, data, yname, tname, idname, gname,
                         xformla, outcome_model, estimand, control_group,
                         doubly_robust, boot_type, anticipation, glist, tlist) {

  if (boot_type == "multiplier") {
    # Multiplier (wild) bootstrap: resample weights
    weights_vec <- stats::rexp(n_units)
    weights_vec <- weights_vec / mean(weights_vec)
    id_weight_map <- stats::setNames(weights_vec, ids)
    boot_data <- data
    boot_data$.boot_wt <- id_weight_map[as.character(boot_data[[idname]])]
  } else {
    # Empirical bootstrap: resample units with replacement
    sampled_ids <- sample(ids, n_units, replace = TRUE)
    boot_data   <- do.call(rbind, lapply(seq_along(sampled_ids), function(k) {
      sub <- data[data[[idname]] == sampled_ids[k], , drop = FALSE]
      sub[[idname]] <- paste0("boot_", k)
      sub
    }))
  }

  atts <- numeric(length(glist) * length(tlist))
  idx  <- 0L

  for (g in glist) {
    for (t in tlist) {
      idx <- idx + 1L
      pre_period <- g - 1L - anticipation
      res <- tryCatch(
        .compute_attgt_single(
          data          = boot_data,
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
        ),
        error = function(e) list(att = NA_real_)
      )
      atts[idx] <- if (is.null(res$att)) NA_real_ else res$att
    }
  }
  atts
}


# ============================================================
# Sandwich / Delta-Method SEs (no bootstrap)
# ============================================================

.sandwich_se <- function(attgt_df, data, yname, tname, idname, gname,
                         xformla, outcome_model, alpha) {

  # Analytical SEs via influence function approach
  # For each (g,t), estimate SE using the empirical influence function
  # of the DR estimator.

  ids     <- unique(data[[idname]])
  n       <- length(ids)

  for (i in seq_len(nrow(attgt_df))) {
    g <- attgt_df$group[i]
    t <- attgt_df$time[i]

    if (is.na(attgt_df$att[i])) next

    # Compute influence function for DR estimator at this (g,t)
    se_hat <- tryCatch({
      .if_se(data = data, yname = yname, tname = tname, idname = idname,
             gname = gname, xformla = xformla, g = g, t = t,
             outcome_model = outcome_model, n = n)
    }, error = function(e) NA_real_)

    attgt_df$se[i] <- se_hat
  }

  attgt_df$ci_lo <- attgt_df$att - stats::qnorm(1 - alpha / 2) * attgt_df$se
  attgt_df$ci_hi <- attgt_df$att + stats::qnorm(1 - alpha / 2) * attgt_df$se
  attgt_df
}


.if_se <- function(data, yname, tname, idname, gname, xformla, g, t,
                   outcome_model, n) {

  pre_period <- g - 1L

  # Get the relevant 2-period subset
  relevant_ids <- data[[idname]][data[[gname]] == g | data[[gname]] == 0]
  sub <- data[data[[idname]] %in% relevant_ids &
                data[[tname]] %in% c(pre_period, t), , drop = FALSE]

  if (nrow(sub) < 4) return(NA_real_)

  sub$D <- as.integer(sub[[gname]] == g)
  sub_wide <- .make_wide(sub, idname, tname, yname, pre_period, t)
  sub_wide$D <- as.integer(sub_wide[[paste0("D_t1")]])

  if (nrow(sub_wide) < 4) return(NA_real_)

  Y0   <- sub_wide[[paste0(yname, "_t0")]]
  Y1   <- sub_wide[[paste0(yname, "_t1")]]
  D    <- sub_wide$D
  n_g  <- nrow(sub_wide)

  # Propensity score
  ps_formula <- stats::update(xformla, D ~ .)
  ps_fit     <- suppressWarnings(stats::glm(ps_formula, data = sub_wide,
                                            family = stats::binomial()))
  pscore     <- pmin(pmax(stats::predict(ps_fit, type = "response"), 1e-6), 1 - 1e-6)

  # Control outcome model
  controls_df   <- sub_wide[D == 0, , drop = FALSE]
  controls_df$DY <- Y1[D == 0] - Y0[D == 0]
  or_formula    <- stats::update(xformla, DY ~ .)
  or_fit        <- suppressWarnings(stats::lm(or_formula, data = controls_df))
  mu_hat        <- stats::predict(or_fit, newdata = sub_wide)

  DeltaY <- Y1 - Y0
  p_D    <- mean(D)
  w_t    <- D / p_D
  w_c    <- (1 - D) * pscore / (1 - pscore) / p_D

  # Influence function of DR-ATT
  psi_att <- w_t * (DeltaY - mu_hat) - w_c * (DeltaY - mu_hat)

  # SE via IF variance
  se_hat <- sqrt(stats::var(psi_att, na.rm = TRUE) / n_g) * sqrt(n_g / n)
  se_hat
}
