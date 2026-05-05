#' @title Aggregate ATT(g,t) Estimates for Nonlinear DiD
#'
#' @description
#' Aggregates the group-time average treatment effects from
#' \code{\link{nonlinear_attgt}} into interpretable summary parameters.
#' Provides event-study (dynamic), group-level, calendar-time, and
#' overall ATT aggregations - each appropriate for nonlinear settings.
#'
#' @param obj An object of class \code{nonlinear_attgt} from
#'   \code{\link{nonlinear_attgt}}.
#' @param type Character. The aggregation type:
#'   \itemize{
#'     \item \code{"dynamic"}: Event-study / dynamic treatment effects.
#'       Averages ATT(g,t) across groups g for each relative time e = t - g.
#'     \item \code{"group"}: Group-specific ATT. Averages over post-treatment
#'       periods within each treated cohort g.
#'     \item \code{"calendar"}: Calendar-time ATT. Averages over groups for
#'       each calendar time t.
#'     \item \code{"simple"}: Overall average ATT, weighted by cohort size.
#'   }
#' @param na.rm Logical. Remove NA ATT(g,t) estimates. Default TRUE.
#' @param min_periods Integer. Minimum number of ATT(g,t) observations
#'   required for an aggregated estimate to be reported. Default 1.
#' @param weights Character. Weighting scheme for aggregation:
#'   \itemize{
#'     \item \code{"equal"}: Equal-weight across (g,t) cells (default).
#'     \item \code{"sample"}: Weight by treated sample size in each (g,t).
#'   }
#'
#' @return An object of class \code{nonlinear_aggte} with slots:
#'   \describe{
#'     \item{agg}{Data frame with aggregated ATT, SE, and CI.}
#'     \item{type}{The aggregation type used.}
#'     \item{overall_att}{Scalar overall ATT estimate.}
#'     \item{overall_se}{SE for overall ATT.}
#'   }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat  <- sim_binary_panel(n = 400, nperiods = 8, prop_treated = 0.5)
#' res  <- nonlinear_attgt(dat, yname = "y", tname = "period",
#'                          idname = "id", gname = "g",
#'                          outcome_model = "logit")
#' agg  <- nonlinear_aggte(res, type = "dynamic")
#' plot(agg)
#'
#' }
#' @export
nonlinear_aggte <- function(
    obj,
    type    = c("dynamic", "group", "calendar", "simple"),
    na.rm   = TRUE,
    min_periods = 1L,
    weights = c("equal", "sample")
) {
  if (!inherits(obj, "nonlinear_attgt")) {
    stop("'obj' must be an object of class 'nonlinear_attgt'.")
  }

  type    <- match.arg(type)
  weights <- match.arg(weights)

  attgt   <- obj$attgt
  if (na.rm) attgt <- attgt[!is.na(attgt$att), ]

  if (nrow(attgt) == 0) stop("No valid ATT(g,t) estimates available.")

  # Compute weights
  attgt$wt <- if (weights == "sample") {
    attgt$n_treated / sum(attgt$n_treated, na.rm = TRUE)
  } else {
    rep(1 / nrow(attgt), nrow(attgt))
  }

  agg_df <- switch(type,
    "dynamic" = {
      attgt$e <- attgt$time - attgt$group
      groups  <- split(attgt, attgt$e)
      do.call(rbind, lapply(names(groups), function(e_val) {
        sub <- groups[[e_val]]
        if (nrow(sub) < min_periods) return(NULL)
        w   <- sub$wt / sum(sub$wt)
        est <- sum(w * sub$att)
        se  <- if (!all(is.na(sub$se))) sqrt(sum(w^2 * sub$se^2, na.rm = TRUE)) else NA_real_
        data.frame(
          label   = as.numeric(e_val),
          att     = est,
          se      = se,
          ci_lo   = est - stats::qnorm(0.975) * se,
          ci_hi   = est + stats::qnorm(0.975) * se,
          n_groups = nrow(sub),
          post    = as.numeric(e_val) >= 0
        )
      }))
    },

    "group" = {
      groups <- split(attgt[attgt$post, ], attgt$group[attgt$post])
      do.call(rbind, lapply(names(groups), function(g_val) {
        sub <- groups[[g_val]]
        if (nrow(sub) < min_periods) return(NULL)
        w   <- sub$wt / sum(sub$wt)
        est <- sum(w * sub$att)
        se  <- if (!all(is.na(sub$se))) sqrt(sum(w^2 * sub$se^2, na.rm = TRUE)) else NA_real_
        data.frame(
          label    = as.numeric(g_val),
          att      = est,
          se       = se,
          ci_lo    = est - stats::qnorm(0.975) * se,
          ci_hi    = est + stats::qnorm(0.975) * se,
          n_periods = nrow(sub),
          post     = TRUE
        )
      }))
    },

    "calendar" = {
      post_df <- attgt[attgt$post, ]
      groups  <- split(post_df, post_df$time)
      do.call(rbind, lapply(names(groups), function(t_val) {
        sub <- groups[[t_val]]
        if (nrow(sub) < min_periods) return(NULL)
        w   <- sub$wt / sum(sub$wt)
        est <- sum(w * sub$att)
        se  <- if (!all(is.na(sub$se))) sqrt(sum(w^2 * sub$se^2, na.rm = TRUE)) else NA_real_
        data.frame(
          label   = as.numeric(t_val),
          att     = est,
          se      = se,
          ci_lo   = est - stats::qnorm(0.975) * se,
          ci_hi   = est + stats::qnorm(0.975) * se,
          n_groups = nrow(sub),
          post     = TRUE
        )
      }))
    },

    "simple" = {
      post_df <- attgt[attgt$post, ]
      w   <- post_df$wt / sum(post_df$wt)
      est <- sum(w * post_df$att)
      se  <- if (!all(is.na(post_df$se))) {
        sqrt(sum(w^2 * post_df$se^2, na.rm = TRUE))
      } else NA_real_
      data.frame(
        label    = "Overall",
        att      = est,
        se       = se,
        ci_lo    = est - stats::qnorm(0.975) * se,
        ci_hi    = est + stats::qnorm(0.975) * se,
        n_groups = nrow(post_df),
        post     = TRUE
      )
    }
  )

  if (is.null(agg_df) || nrow(agg_df) == 0) {
    stop("Aggregation produced no results. Try adjusting 'min_periods'.")
  }

  # Overall ATT
  if (type == "simple") {
    overall_att <- agg_df$att[1]
    overall_se  <- agg_df$se[1]
  } else {
    post_df <- attgt[attgt$post, ]
    w   <- post_df$wt / sum(post_df$wt)
    overall_att <- sum(w * post_df$att)
    overall_se  <- if (!all(is.na(post_df$se))) {
      sqrt(sum(w^2 * post_df$se^2, na.rm = TRUE))
    } else NA_real_
  }

  out <- list(
    agg         = agg_df,
    type        = type,
    overall_att = overall_att,
    overall_se  = overall_se,
    args        = obj$args
  )
  class(out) <- "nonlinear_aggte"
  out
}
