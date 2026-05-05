#' @title S3 Methods for NonlinearDiD Objects
#' @description Print, summary, and plot methods for \code{nonlinear_attgt}
#'   and \code{nonlinear_aggte} objects.
#' @name nonlineardid_methods
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline
#'   geom_vline facet_wrap theme_bw theme element_text labs
#'   scale_color_manual geom_ribbon geom_line
NULL


# ============================================================
# print methods
# ============================================================

#' @export
print.nonlinear_attgt <- function(x, ...) {
  cat("\n=== Nonlinear Staggered DiD: Group-Time ATT(g,t) ===\n\n")
  cat("Outcome model:  ", x$args$outcome_model, "\n")
  cat("Estimand:       ", x$args$estimand, "\n")
  cat("Control group:  ", x$args$control_group, "\n")
  cat("Groups (cohorts):", paste(x$args$glist, collapse = ", "), "\n")
  cat("Time periods:   ", paste(x$args$tlist, collapse = ", "), "\n\n")

  att <- x$attgt
  n_valid    <- sum(!is.na(att$att))
  n_post     <- sum(att$post & !is.na(att$att))
  n_pre      <- n_valid - n_post
  n_converge <- sum(att$converged, na.rm = TRUE)

  cat(sprintf("ATT(g,t) estimates: %d total (%d pre, %d post), %d converged\n\n",
              n_valid, n_pre, n_post, n_converge))

  # Show first few rows
  cols_show <- intersect(c("group", "time", "att", "se", "ci_lo", "ci_hi", "post"),
                         names(att))
  print(utils::head(att[, cols_show, drop = FALSE], 12), digits = 4, row.names = FALSE)

  if (nrow(att) > 12) cat("... (", nrow(att) - 12, "more rows)\n")
  cat("\nUse summary() for overall ATT. Use nonlinear_aggte() for aggregation.\n")
  invisible(x)
}


#' @export
print.nonlinear_aggte <- function(x, ...) {
  cat("\n=== Nonlinear DiD: Aggregated ATT ===\n\n")
  cat("Aggregation type:", x$type, "\n")
  cat("Outcome model:  ", x$args$outcome_model, "\n\n")

  if (!is.na(x$overall_att)) {
    z     <- x$overall_att / x$overall_se
    pval  <- 2 * stats::pnorm(-abs(z))
    cat(sprintf("Overall ATT: %.4f  (SE: %.4f, t = %.3f, p = %.4f)\n",
                x$overall_att, x$overall_se, z, pval))
    cat(sprintf("95%% CI: [%.4f, %.4f]\n\n",
                x$overall_att - 1.96 * x$overall_se,
                x$overall_att + 1.96 * x$overall_se))
  }

  print(x$agg, digits = 4, row.names = FALSE)
  invisible(x)
}


# ============================================================
# summary methods
# ============================================================

#' @export
summary.nonlinear_attgt <- function(object, ...) {
  att <- object$attgt

  post_df <- att[att$post & !is.na(att$att), ]
  pre_df  <- att[!att$post & !is.na(att$att), ]

  cat("\n=== Nonlinear DiD Summary ===\n\n")
  cat("Model:        ", object$args$outcome_model, "\n")
  cat("Estimand:     ", object$args$estimand, "\n")
  cat("Control group:", object$args$control_group, "\n\n")

  cat("--- Post-treatment ATT(g,t) ---\n")
  if (nrow(post_df) > 0) {
    cat(sprintf("  Mean ATT: %.4f\n", mean(post_df$att, na.rm = TRUE)))
    cat(sprintf("  Median ATT: %.4f\n", stats::median(post_df$att, na.rm = TRUE)))
    cat(sprintf("  Min ATT: %.4f  |  Max ATT: %.4f\n",
                min(post_df$att), max(post_df$att)))
    cat(sprintf("  N(g,t) cells: %d\n", nrow(post_df)))
  } else cat("  None.\n")

  cat("\n--- Pre-treatment ATT(g,t) (should be ~0 under PT) ---\n")
  if (nrow(pre_df) > 0) {
    cat(sprintf("  Mean: %.4f  |  SD: %.4f\n",
                mean(pre_df$att, na.rm = TRUE), stats::sd(pre_df$att, na.rm = TRUE)))
    cat(sprintf("  N(g,t) cells: %d\n", nrow(pre_df)))
  } else cat("  None.\n")

  cat("\n--- Convergence ---\n")
  cat(sprintf("  Converged: %d / %d\n",
              sum(att$converged, na.rm = TRUE), nrow(att)))

  invisible(object)
}


#' @export
summary.nonlinear_aggte <- function(object, ...) {
  print(object)
}


# ============================================================
# plot methods
# ============================================================

#' @title Plot ATT(g,t) Estimates
#'
#' @description Produces a faceted scatter plot of ATT(g,t) estimates with
#'   confidence intervals, one panel per treatment cohort.
#'
#' @param x An object of class \code{nonlinear_attgt}.
#' @param ... Additional arguments (unused).
#' @param alpha Numeric. Significance level for CI. Default 0.05.
#' @param point_size Numeric. Size of estimate points. Default 2.
#'
#' @return A \code{ggplot2} object.
#' @export
plot.nonlinear_attgt <- function(x, ..., alpha = 0.05, point_size = 2) {
  att <- x$attgt
  att <- att[!is.na(att$att), ]

  if (!all(c("ci_lo", "ci_hi") %in% names(att))) {
    att$ci_lo <- att$att - 1.96 * att$se
    att$ci_hi <- att$att + 1.96 * att$se
  }

  att$cohort_label <- paste0("Cohort g = ", att$group)
  att$period_type  <- ifelse(att$post, "Post-treatment", "Pre-treatment")

  p <- ggplot2::ggplot(att, ggplot2::aes(x = time, y = att, color = period_type)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = group - 0.5), linetype = "dotted",
               color = "grey40", linewidth = 0.5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lo, ymax = ci_hi), width = 0.2,
                  alpha = 0.7) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::facet_wrap(~ cohort_label, scales = "free_x") +
    ggplot2::scale_color_manual(
      values = c("Post-treatment" = "#2166AC", "Pre-treatment" = "#D73027"),
      name   = ""
    ) +
    ggplot2::labs(
      title    = paste0("ATT(g,t) Estimates - Nonlinear DiD (",
                        x$args$outcome_model, ")"),
      subtitle = paste0("Estimand: ", x$args$estimand,
                        " | Control: ", x$args$control_group),
      x        = "Calendar Time",
      y        = paste0("ATT [", x$args$estimand, "]"),
      caption  = "Dashed: zero line. Dotted: first treatment period."
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position  = "bottom",
      strip.text       = ggplot2::element_text(size = 10, face = "bold"),
      plot.title       = ggplot2::element_text(size = 12, face = "bold"),
      axis.text        = ggplot2::element_text(size = 9)
    )

  p
}


#' @title Plot Aggregated DiD Estimates
#'
#' @description Plots event-study, group-level, calendar, or overall
#'   aggregated ATT estimates from \code{\link{nonlinear_aggte}}.
#'
#' @param x An object of class \code{nonlinear_aggte}.
#' @param ... Additional arguments (unused).
#'
#' @return A \code{ggplot2} object.
#' @export
plot.nonlinear_aggte <- function(x, ...) {
  agg  <- x$agg
  type <- x$type

  if (!("ci_lo" %in% names(agg))) {
    agg$ci_lo <- agg$att - 1.96 * agg$se
    agg$ci_hi <- agg$att + 1.96 * agg$se
  }

  title_map <- c(
    dynamic  = "Event-Study (Dynamic) ATT",
    group    = "Group-Specific ATT",
    calendar = "Calendar-Time ATT",
    simple   = "Overall ATT"
  )

  xlab_map  <- c(
    dynamic  = "Relative Time to Treatment",
    group    = "Treatment Cohort (g)",
    calendar = "Calendar Time",
    simple   = ""
  )

  if (type == "dynamic") {
    agg$post  <- agg$label >= 0
    agg$color <- ifelse(agg$post, "#2166AC", "#D73027")

    p <- ggplot2::ggplot(agg, ggplot2::aes(x = label, y = att)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_vline(xintercept = -0.5, linetype = "dotted", color = "grey40") +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.15,
                  fill = "#2166AC") +
      ggplot2::geom_line(color = "#2166AC", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(color = post), size = 2.5) +
      ggplot2::scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D73027"),
                         labels = c("TRUE" = "Post", "FALSE" = "Pre"),
                         name = "") +
      ggplot2::labs(
        title    = title_map[type],
        subtitle = paste0("Model: ", x$args$outcome_model,
                          " | Overall ATT = ",
                          round(x$overall_att, 4)),
        x        = xlab_map[type],
        y        = "ATT",
        caption  = "Shaded region: 95% pointwise CI. Dotted: event time -0.5."
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom",
            plot.title = ggplot2::element_text(face = "bold"))
  } else {
    p <- ggplot2::ggplot(agg, ggplot2::aes(x = factor(label), y = att)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lo, ymax = ci_hi), width = 0.25,
                    color = "#2166AC") +
      ggplot2::geom_point(size = 3, color = "#2166AC") +
      ggplot2::labs(
        title = title_map[type],
        subtitle = paste0("Model: ", x$args$outcome_model),
        x     = xlab_map[type],
        y     = "ATT"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
  }

  p
}


# ============================================================
# print methods for 2x2 estimators
# ============================================================

#' @export
print.binary_did_logit <- function(x, ...) {
  .print_2x2_did(x, "Binary DiD (Logit)", "log-odds")
}

#' @export
print.binary_did_probit <- function(x, ...) {
  .print_2x2_did(x, "Binary DiD (Probit)", "probit index")
}

#' @export
print.binary_did_dr <- function(x, ...) {
  cat("\n=== Doubly-Robust Binary DiD ===\n\n")
  cat(sprintf("ATT:   %.4f  (SE: %.4f)\n", x$att, x$se))
  cat(sprintf("95%% CI: [%.4f, %.4f]\n", x$ci_lo, x$ci_hi))
  cat(sprintf("t-stat: %.3f  |  p-value: %.4f\n", x$tstat, x$pval))
  cat(sprintf("N (treated): %d  |  N (control): %d\n", x$n_treat, x$n_cont))
  invisible(x)
}

#' @export
print.count_did_poisson <- function(x, ...) {
  cat("\n=== Count DiD (Poisson QMLE) ===\n\n")
  cat(sprintf("ATT (log rate ratio): %.4f  (SE: %.4f)\n", x$att_log_rr, x$se_log_rr))
  cat(sprintf("Rate ratio:           %.4f  [95%% CI: %.4f, %.4f]\n",
              x$rate_ratio, x$ci_lo_rr, x$ci_hi_rr))
  cat(sprintf("ATT (APE, count scale): %.4f\n", x$att_ape))
  cat(sprintf("t-stat: %.3f  |  p-value: %.4f\n", x$tstat, x$pval))
  invisible(x)
}

#' @export
print.odds_ratio_did <- function(x, ...) {
  cat("\n=== Odds-Ratio DiD ===\n\n")
  cat("Interpretation:", x$interpretation, "\n\n")
  cat(sprintf("Log(OR-DiD): %.4f  (SE: %.4f)\n", x$log_or_did, x$se_log_or))
  cat(sprintf("95%% CI [OR]: [%.4f, %.4f]\n", x$ci_lo_or, x$ci_hi_or))
  cat(sprintf("t-stat: %.3f  |  p-value: %.4f\n", x$tstat, x$pval))
  invisible(x)
}


.print_2x2_did <- function(x, title, scale_name) {
  cat(paste0("\n=== ", title, " ===\n\n"))
  cat(sprintf("DiD on %s scale:  %.4f  (SE: %.4f)\n",
              scale_name, x$att_link, x$se_link))
  cat(sprintf("95%% CI (%s):    [%.4f, %.4f]\n",
              scale_name, x$ci_lo_link, x$ci_hi_link))
  cat(sprintf("Average partial effect (prob): %.4f\n", x$att_ape))
  cat(sprintf("t-stat: %.3f  |  p-value: %.4f\n", x$tstat, x$pval))
  invisible(x)
}

# Suppress R CMD check NOTEs for ggplot2 NSE column names
utils::globalVariables(c(
  "label", "att", "ci_lo", "ci_hi", "post",
  "time", "period_type", "group",
  "cohort_label"
))
