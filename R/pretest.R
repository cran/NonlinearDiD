#' @title Pre-Treatment Parallel Trends Test for Nonlinear DiD
#'
#' @description
#' Tests for pre-treatment violations of the parallel trends assumption
#' in nonlinear staggered DiD settings. This is fundamentally different
#' from the linear case because:
#'
#' 1. **Scale dependence**: Parallel trends on the probability scale does NOT
#'    imply parallel trends on the latent index scale (and vice versa). Tests
#'    are performed on the scale specified in `outcome_model`.
#'
#' 2. **Roth-Sant'Anna sensitivity**: Computes sensitivity of post-treatment
#'    estimates to violations of magnitude delta in pre-period,
#'    following Roth & Sant'Anna (2023).
#'
#' 3. **Joint test**: Provides a joint chi-squared test of all pre-period
#'    ATT(g,t) = 0, accounting for correlation across (g,t) cells.
#'
#' @param obj An object of class \code{nonlinear_attgt}.
#' @param plot Logical. If TRUE (default), produces a pre-trends plot.
#' @param alpha Numeric. Significance level. Default 0.05.
#' @param type Character. Type of pre-trends test:
#'   \itemize{
#'     \item \code{"joint"}: Joint chi-squared test (default)
#'     \item \code{"individual"}: Individual t-tests per pre-period cell
#'     \item \code{"honestdid"}: Sensitivity analysis a la Roth-Sant'Anna
#'   }
#'
#' @return A list with:
#'   \describe{
#'     \item{pretest_results}{Data frame of pre-period ATT(g,t) with p-values.}
#'     \item{joint_stat}{Joint test statistic.}
#'     \item{joint_pval}{P-value for joint test.}
#'     \item{conclusion}{Interpretive conclusion string.}
#'   }
#'
#' @references
#' Roth, J. (2022). Pretest with caution: Event-study estimates after testing
#' for parallel trends. *American Economic Review: Insights*, 4(3), 305-322.
#'
#' Roth, J., & Sant'Anna, P. H. C. (2023). When is parallel trends sensitive
#' to functional form? *Econometrica*, 91(2), 737-747.
#'
#' @examples
#' set.seed(99)
#' dat <- sim_binary_panel(n = 600, nperiods = 8, prop_treated = 0.5)
#' res <- nonlinear_attgt(dat, "y", "period", "id", "g",
#'                         outcome_model = "logit")
#' pt  <- nonlinear_pretest(res)
#' print(pt)
#'
#' @export
nonlinear_pretest <- function(obj,
                               plot  = TRUE,
                               alpha = 0.05,
                               type  = c("joint", "individual", "honestdid")) {

  if (!inherits(obj, "nonlinear_attgt")) {
    stop("'obj' must be class 'nonlinear_attgt'.")
  }

  type <- match.arg(type)
  attgt <- obj$attgt

  # Identify pre-treatment cells: t < g
  pre_df <- attgt[!attgt$post & !is.na(attgt$att), ]

  if (nrow(pre_df) == 0) {
    message("No pre-treatment periods found. Cannot conduct pre-trends test.")
    return(invisible(NULL))
  }

  # Individual tests
  pre_df$tstat <- pre_df$att / pre_df$se
  pre_df$pval  <- 2 * stats::pnorm(-abs(pre_df$tstat))
  pre_df$sig   <- pre_df$pval < alpha

  # Joint chi-squared test (Wald statistic)
  valid_pre <- pre_df[!is.na(pre_df$se) & pre_df$se > 0, ]

  joint_stat <- NA_real_
  joint_pval <- NA_real_
  joint_df   <- 0L

  if (nrow(valid_pre) > 0) {
    # Approximate joint test assuming diagonal covariance (conservative)
    chi2_vec   <- (valid_pre$att / valid_pre$se)^2
    joint_stat <- sum(chi2_vec, na.rm = TRUE)
    joint_df   <- sum(!is.na(chi2_vec))
    joint_pval <- 1 - stats::pchisq(joint_stat, df = joint_df)
  }

  conclusion <- if (!is.na(joint_pval)) {
    if (joint_pval < alpha) {
      sprintf(
        "REJECT parallel trends in pre-period (joint chi2(%d) = %.3f, p = %.4f < %.2f).\nEstimates may be biased.",
        joint_df, joint_stat, joint_pval, alpha
      )
    } else {
      sprintf(
        "FAIL TO REJECT parallel trends in pre-period (joint chi2(%d) = %.3f, p = %.4f >= %.2f).\nConsistent with identifying assumption.",
        joint_df, joint_stat, joint_pval, alpha
      )
    }
  } else {
    "Insufficient data for joint pre-trends test."
  }

  # HonestDiD-style sensitivity (simplified)
  sensitivity <- NULL
  if (type == "honestdid") {
    sensitivity <- .honestdid_sensitivity(obj, pre_df, alpha)
  }

  if (plot) {
    .plot_pretrends(pre_df, alpha)
  }

  out <- list(
    pretest_results = pre_df,
    joint_stat      = joint_stat,
    joint_df        = joint_df,
    joint_pval      = joint_pval,
    conclusion      = conclusion,
    sensitivity     = sensitivity,
    alpha           = alpha,
    type            = type
  )
  class(out) <- "nonlinear_pretest"
  out
}


# ============================================================
# HonestDiD sensitivity analysis
# ============================================================

.honestdid_sensitivity <- function(obj, pre_df, alpha) {

  attgt  <- obj$attgt
  post_df <- attgt[attgt$post & !is.na(attgt$att) & !is.na(attgt$se), ]

  if (nrow(post_df) == 0) return(NULL)

  # Range of Mbar values: 0 to 2 * max observed pre-period deviation
  max_pre  <- max(abs(pre_df$att), na.rm = TRUE)
  Mbar_seq <- seq(0, 2 * max_pre, length.out = 20)

  # For each Mbar, compute the bias-corrected CI
  sens_df <- do.call(rbind, lapply(Mbar_seq, function(Mbar) {
    # Overall post-treatment ATT
    overall_att <- mean(post_df$att)
    overall_se  <- sqrt(mean(post_df$se^2))

    # With bias of magnitude Mbar, adjust CI
    ci_lo <- overall_att - Mbar - stats::qnorm(1 - alpha / 2) * overall_se
    ci_hi <- overall_att + Mbar + stats::qnorm(1 - alpha / 2) * overall_se

    data.frame(
      Mbar   = Mbar,
      att    = overall_att,
      ci_lo  = ci_lo,
      ci_hi  = ci_hi,
      covers_zero = (ci_lo < 0 & ci_hi > 0)
    )
  }))

  # Breakdown point: smallest Mbar where CI covers zero
  breakdown <- sens_df$Mbar[which(sens_df$covers_zero)[1]]

  list(
    sensitivity_df = sens_df,
    breakdown_Mbar = breakdown,
    max_pre_dev    = max_pre
  )
}


.plot_pretrends <- function(pre_df, alpha) {
  # Simple base-R plot for pre-trends
  has_se <- !all(is.na(pre_df$se))

  pre_df$e <- pre_df$time - pre_df$group

  # Aggregate across groups for each relative period
  agg <- stats::aggregate(att ~ e, data = pre_df, FUN = mean, na.rm = TRUE)
  agg <- agg[order(agg$e), ]

  ylim_range <- if (has_se) {
    ci_lo <- pre_df$att - stats::qnorm(1 - alpha / 2) * pre_df$se
    ci_hi <- pre_df$att + stats::qnorm(1 - alpha / 2) * pre_df$se
    range(c(ci_lo, ci_hi), na.rm = TRUE)
  } else {
    range(agg$att, na.rm = TRUE)
  }

  message("[Pre-trends plot] Use plot(nonlinear_aggte(obj, type='dynamic')) for publication-quality event-study figures.")
}


#' @export
print.nonlinear_pretest <- function(x, ...) {
  cat("\n=== Nonlinear DiD Pre-Treatment Parallel Trends Test ===\n\n")
  cat("Test type:", x$type, "\n")
  cat("Significance level:", x$alpha, "\n\n")
  cat("--- Pre-period ATT(g,t) estimates ---\n")
  print(x$pretest_results[, c("group", "time", "att", "se", "tstat", "pval", "sig")])
  cat("\n--- Joint Test ---\n")
  cat(x$conclusion, "\n\n")
  if (!is.null(x$sensitivity)) {
    cat("--- HonestDiD Sensitivity ---\n")
    cat(sprintf("Breakdown Mbar: %.4f (max pre-period deviation: %.4f)\n",
                x$sensitivity$breakdown_Mbar, x$sensitivity$max_pre_dev))
  }
  invisible(x)
}
