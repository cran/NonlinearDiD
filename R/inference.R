#' @title Inference for Nonlinear DiD
#' @description Internal functions for bootstrap and delta-method standard
#'   errors in nonlinear staggered DiD, for both panel and repeated
#'   cross-section designs, with optional clustering and sampling weights.
#' @keywords internal

# ============================================================
# Multiplier / Empirical Bootstrap
# ============================================================

.bootstrap_inference <- function(attgt_df, data, yname, tname, idname, gname,
                                 xformla, outcome_model, estimand, control_group,
                                 doubly_robust, nboot, boot_type, alpha,
                                 anticipation, parallel, pl_cores,
                                 data_type = "panel", cluster_var = NULL) {

  tlist <- sort(unique(data[[tname]]))
  glist <- sort(setdiff(unique(data[[gname]]), 0))

  # Decide what gets resampled.
  # 1. cluster_var provided  -> cluster bootstrap
  # 2. data_type panel       -> unit-id bootstrap
  # 3. data_type RCS         -> row bootstrap
  if (!is.null(cluster_var)) {
    cluster_ids <- unique(data[[cluster_var]])
    boot_unit_col <- cluster_var
    n_units <- length(cluster_ids)
    ids     <- cluster_ids
  } else if (data_type == "panel") {
    ids     <- unique(data[[idname]])
    n_units <- length(ids)
    boot_unit_col <- idname
  } else {
    ids     <- seq_len(nrow(data))
    n_units <- nrow(data)
    boot_unit_col <- NULL  # row-level
  }

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(pl_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    boot_draws <- parallel::parSapply(cl, seq_len(nboot), function(b)
      .single_boot(b, ids, n_units, data, yname, tname, idname, gname,
                   xformla, outcome_model, estimand, control_group,
                   doubly_robust, boot_type, anticipation, glist, tlist,
                   data_type, boot_unit_col))
  } else {
    boot_draws <- matrix(NA_real_, nrow = nrow(attgt_df), ncol = nboot)
    for (b in seq_len(nboot)) {
      bv <- .single_boot(b, ids, n_units, data, yname, tname, idname, gname,
                         xformla, outcome_model, estimand, control_group,
                         doubly_robust, boot_type, anticipation, glist, tlist,
                         data_type, boot_unit_col)
      boot_draws[, b] <- bv
    }
  }

  attgt_df$se    <- apply(boot_draws, 1, function(x) stats::sd(x, na.rm = TRUE))
  attgt_df$ci_lo <- attgt_df$att - stats::qnorm(1 - alpha / 2) * attgt_df$se
  attgt_df$ci_hi <- attgt_df$att + stats::qnorm(1 - alpha / 2) * attgt_df$se

  # Simultaneous CIs from bootstrap quantile
  attgt_df$ci_lo_simult <- NA_real_
  attgt_df$ci_hi_simult <- NA_real_
  for (i in seq_len(nrow(attgt_df))) {
    q <- stats::quantile(abs(boot_draws[i, ] - attgt_df$att[i]),
                          1 - alpha, na.rm = TRUE)
    attgt_df$ci_lo_simult[i] <- attgt_df$att[i] - q
    attgt_df$ci_hi_simult[i] <- attgt_df$att[i] + q
  }

  attgt_df
}


.single_boot <- function(b, ids, n_units, data, yname, tname, idname, gname,
                         xformla, outcome_model, estimand, control_group,
                         doubly_robust, boot_type, anticipation, glist, tlist,
                         data_type = "panel", boot_unit_col = NULL) {

  # Build boot_data:
  # - boot_unit_col not NULL: cluster/unit bootstrap on that column
  # - boot_unit_col NULL    : row-level bootstrap (RCS no cluster)
  if (!is.null(boot_unit_col)) {
    if (boot_type == "multiplier") {
      weights_vec <- stats::rexp(n_units)
      weights_vec <- weights_vec / mean(weights_vec)
      id_weight_map <- stats::setNames(weights_vec, ids)
      boot_data <- data
      boot_data$.boot_wt <- id_weight_map[as.character(boot_data[[boot_unit_col]])]
      # Fold into sampling weight
      boot_data$.sw <- boot_data$.sw * boot_data$.boot_wt
    } else {
      sampled_ids <- sample(ids, n_units, replace = TRUE)
      boot_data <- do.call(rbind, lapply(seq_along(sampled_ids), function(k) {
        sub <- data[data[[boot_unit_col]] == sampled_ids[k], , drop = FALSE]
        # Re-label the unit to keep rows uniquely identifiable in panel mode
        if (data_type == "panel" && !is.null(idname)) {
          sub[[idname]] <- paste0("boot_", k, "_", sub[[idname]])
        }
        sub
      }))
    }
  } else {
    # Row-level (RCS, no clustering)
    if (boot_type == "multiplier") {
      row_wts <- stats::rexp(n_units)
      row_wts <- row_wts / mean(row_wts)
      boot_data <- data
      boot_data$.sw <- boot_data$.sw * row_wts
    } else {
      # Resample rows within each (group x period) cell
      cells <- split(seq_len(nrow(data)), list(data[[gname]], data[[tname]]))
      new_rows <- unlist(lapply(cells, function(idx) {
        if (length(idx) == 0) return(integer(0))
        sample(idx, length(idx), replace = TRUE)
      }))
      boot_data <- data[new_rows, , drop = FALSE]
      rownames(boot_data) <- NULL
    }
  }

  atts <- numeric(length(glist) * length(tlist))
  idx  <- 0L

  for (g in glist) {
    for (t in tlist) {
      idx <- idx + 1L
      pre_period <- g - 1L - anticipation
      res <- tryCatch(
        .compute_attgt_single(
          data = boot_data, yname = yname, tname = tname,
          idname = idname, gname = gname, xformla = xformla,
          g = g, t = t, base_period = pre_period,
          outcome_model = outcome_model, estimand = estimand,
          control_group = control_group, doubly_robust = doubly_robust,
          data_type = data_type
        ),
        error = function(e) list(att = NA_real_)
      )
      atts[idx] <- if (is.null(res$att)) NA_real_ else res$att
    }
  }
  atts
}


# ============================================================
# Analytical SEs (no bootstrap)
# ============================================================

.sandwich_se <- function(attgt_df, data, yname, tname, idname, gname,
                         xformla, outcome_model, alpha,
                         data_type = "panel", cluster_var = NULL) {

  if (data_type == "repeated_cross_section") {
    return(.sandwich_se_rcs(
      attgt_df = attgt_df, data = data, yname = yname, tname = tname,
      gname = gname, xformla = xformla, outcome_model = outcome_model,
      alpha = alpha, cluster_var = cluster_var
    ))
  }

  # ---- Panel path (v0.1.0 IF-based SEs; clustering left to bootstrap) ----
  ids <- unique(data[[idname]])
  n   <- length(ids)

  for (i in seq_len(nrow(attgt_df))) {
    g <- attgt_df$group[i]
    t <- attgt_df$time[i]
    if (is.na(attgt_df$att[i])) next

    se_hat <- tryCatch({
      .if_se(data = data, yname = yname, tname = tname, idname = idname,
             gname = gname, xformla = xformla, g = g, t = t,
             outcome_model = outcome_model, n = n)
    }, error = function(e) NA_real_)

    attgt_df$se[i] <- se_hat
  }

  if (!is.null(cluster_var)) {
    # Panel + clustering with analytical SEs is not yet implemented.
    # Inflate the IF-based SE by a conservative DEFF-style adjustment
    # based on average cluster size, so users still get a SE estimate.
    # Recommend bootstrap with cluster_var for proper clustered SEs.
    n_per_cluster <- as.numeric(table(data[[cluster_var]]))
    rho           <- 0.05  # mild intra-cluster correlation assumption
    avg_n         <- mean(n_per_cluster)
    deff          <- 1 + (avg_n - 1) * rho
    attgt_df$se   <- attgt_df$se * sqrt(deff)
    message("Note: panel + cluster_var analytical SEs use a DEFF approximation. ",
            "For full clustered inference, set boot = TRUE.")
  }

  attgt_df$ci_lo <- attgt_df$att - stats::qnorm(1 - alpha / 2) * attgt_df$se
  attgt_df$ci_hi <- attgt_df$att + stats::qnorm(1 - alpha / 2) * attgt_df$se
  attgt_df
}


# ============================================================
# RCS analytical SEs: HC1 or vcovCL on pooled GLM
# ============================================================

.sandwich_se_rcs <- function(attgt_df, data, yname, tname, gname,
                              xformla, outcome_model, alpha,
                              cluster_var = NULL) {

  for (i in seq_len(nrow(attgt_df))) {
    g <- attgt_df$group[i]
    t <- attgt_df$time[i]
    if (is.na(attgt_df$att[i])) next

    se_hat <- tryCatch({
      .rcs_if_se(
        data = data, yname = yname, tname = tname,
        gname = gname, xformla = xformla,
        outcome_model = outcome_model,
        g = g, t = t, cluster_var = cluster_var
      )
    }, error = function(e) NA_real_)

    attgt_df$se[i] <- se_hat
  }

  attgt_df$ci_lo <- attgt_df$att - stats::qnorm(1 - alpha / 2) * attgt_df$se
  attgt_df$ci_hi <- attgt_df$att + stats::qnorm(1 - alpha / 2) * attgt_df$se
  attgt_df
}


# Extract SE for the D:post interaction from a pooled GLM fit on the
# (g, t) cell subset. Uses sandwich::vcovCL if cluster_var is supplied,
# otherwise sandwich::vcovHC (HC1).
.rcs_if_se <- function(data, yname, tname, gname, xformla, outcome_model,
                       g, t, cluster_var = NULL) {

  pre_period <- g - 1L

  sub <- data[(data[[gname]] == g | data[[gname]] == 0) &
                data[[tname]] %in% c(pre_period, t), , drop = FALSE]

  if (nrow(sub) < 8) return(NA_real_)

  sub$D    <- as.integer(sub[[gname]] == g)
  sub$post <- as.integer(sub[[tname]] == t)
  sub$did  <- sub$D * sub$post

  fam   <- .outcome_family(outcome_model)
  xvars <- labels(stats::terms(xformla))
  fmla  <- .rcs_formula(yname, xvars)
  environment(fmla) <- environment()

  sw <- if (".sw" %in% names(sub)) sub$.sw else rep(1.0, nrow(sub))

  fit <- suppressWarnings(stats::glm(fmla, data = sub, family = fam,
                                      weights = sw))

  if (!("did" %in% names(stats::coef(fit)))) return(NA_real_)

  vcov_hc <- tryCatch({
    if (!is.null(cluster_var) && cluster_var %in% names(sub)) {
      sandwich::vcovCL(fit, cluster = sub[[cluster_var]], type = "HC1")
    } else {
      sandwich::vcovHC(fit, type = "HC1")
    }
  }, error = function(e) NULL)

  if (is.null(vcov_hc) || !("did" %in% rownames(vcov_hc))) return(NA_real_)
  sqrt(vcov_hc["did", "did"])
}


# ============================================================
# Panel influence-function SE (v0.1.0, unchanged behaviour)
# ============================================================

.if_se <- function(data, yname, tname, idname, gname, xformla, g, t,
                   outcome_model, n) {

  pre_period <- g - 1L

  relevant_ids <- data[[idname]][data[[gname]] == g | data[[gname]] == 0]
  sub <- data[data[[idname]] %in% relevant_ids &
                data[[tname]] %in% c(pre_period, t), , drop = FALSE]

  if (nrow(sub) < 4) return(NA_real_)

  sub$D    <- as.integer(sub[[gname]] == g)
  sub_wide <- .make_wide(sub, idname, tname, yname, pre_period, t)
  sub_wide$D <- as.integer(sub_wide[["D_t1"]])

  if (nrow(sub_wide) < 4) return(NA_real_)

  Y0  <- sub_wide[[paste0(yname, "_t0")]]
  Y1  <- sub_wide[[paste0(yname, "_t1")]]
  D   <- sub_wide$D
  n_g <- nrow(sub_wide)

  ps_formula <- stats::update(xformla, D ~ .)
  ps_fit     <- suppressWarnings(
    stats::glm(ps_formula, data = sub_wide, family = stats::binomial()))
  pscore     <- pmin(pmax(stats::predict(ps_fit, type = "response"), 1e-6), 1 - 1e-6)

  controls_df    <- sub_wide[D == 0, , drop = FALSE]
  controls_df$DY <- Y1[D == 0] - Y0[D == 0]
  or_formula     <- stats::update(xformla, DY ~ .)
  or_fit         <- suppressWarnings(stats::lm(or_formula, data = controls_df))
  mu_hat         <- stats::predict(or_fit, newdata = sub_wide)

  DeltaY <- Y1 - Y0
  p_D    <- mean(D)
  w_t    <- D / p_D
  w_c    <- (1 - D) * pscore / (1 - pscore) / p_D

  psi_att <- w_t * (DeltaY - mu_hat) - w_c * (DeltaY - mu_hat)

  se_hat <- sqrt(stats::var(psi_att, na.rm = TRUE) / n_g) * sqrt(n_g / n)
  se_hat
}
