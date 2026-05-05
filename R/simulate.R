#' @title Simulate Binary Panel Data with Staggered Treatment
#'
#' @description
#' Generates a simulated panel dataset with staggered treatment adoption and
#' a binary outcome. Useful for testing and illustrating nonlinear DiD methods.
#'
#' The data-generating process is:
#' \deqn{Y_{it} = \mathbf{1}\{ \alpha_i + \lambda_t + \delta_{it} \cdot D_{it} + \epsilon_{it} > 0 \}}
#'
#' where \eqn{\alpha_i} is a unit fixed effect, \eqn{\lambda_t} is a time
#' fixed effect, \eqn{\delta_{it}} is the treatment effect (heterogeneous
#' across cohorts), and \eqn{\epsilon_{it}} is logistic noise.
#'
#' @param n Integer. Number of units. Default 500.
#' @param nperiods Integer. Number of time periods. Default 6.
#' @param prop_treated Numeric. Proportion of units ever treated. Default 0.5.
#' @param n_cohorts Integer. Number of treatment cohorts (groups). Default 3.
#' @param true_att Numeric or vector. True ATT for each cohort. Default 0.3.
#' @param base_prob Numeric. Baseline probability P(Y=1) for untreated.
#'   Default 0.3.
#' @param unit_fe_sd Numeric. Std. dev. of unit fixed effects. Default 0.5.
#' @param add_covariates Logical. Add pre-treatment covariates. Default TRUE.
#' @param seed Integer. Random seed. Default NULL.
#'
#' @return A data frame in long format. Columns: \code{id} (unit identifier),
#'   \code{period} (time period 1 to nperiods), \code{y} (binary outcome 0/1),
#'   \code{g} (treatment cohort; 0 = never treated), \code{D} (treatment
#'   indicator), \code{x1} and \code{x2} (covariates, if
#'   \code{add_covariates = TRUE}), and \code{alpha_i} (true unit fixed effect,
#'   for validation).
#'
#' @examples
#' dat <- sim_binary_panel(n = 1000, nperiods = 8, prop_treated = 0.6,
#'                          n_cohorts = 4, true_att = c(0.2, 0.4, 0.3, 0.5))
#' head(dat)
#' table(dat$g)
#'
#' @export
sim_binary_panel <- function(
    n             = 500L,
    nperiods      = 6L,
    prop_treated  = 0.5,
    n_cohorts     = 3L,
    true_att      = 0.3,
    base_prob     = 0.3,
    unit_fe_sd    = 0.5,
    add_covariates = TRUE,
    seed          = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Extend true_att if scalar
  if (length(true_att) == 1) true_att <- rep(true_att, n_cohorts)
  if (length(true_att) < n_cohorts) {
    true_att <- rep(true_att, length.out = n_cohorts)
  }

  n_treat  <- round(n * prop_treated)
  n_never  <- n - n_treat
  n_per_cohort <- floor(n_treat / n_cohorts)

  # Treatment cohorts: evenly spaced after period 1
  cohort_periods <- seq(2, nperiods - 1, length.out = n_cohorts)
  cohort_periods <- round(cohort_periods)
  cohort_periods <- pmin(pmax(cohort_periods, 2), nperiods)

  # Assign units to groups
  group_vec  <- c(rep(0, n_never),
                  unlist(lapply(seq_len(n_cohorts), function(k) {
                    rep(cohort_periods[k], n_per_cohort)
                  })))

  # Pad to n
  while (length(group_vec) < n) group_vec <- c(group_vec, 0L)
  group_vec <- group_vec[seq_len(n)]

  # Cohort ATT map
  cohort_att_map <- stats::setNames(true_att, cohort_periods)

  # Unit fixed effects
  alpha_i <- stats::rnorm(n, mean = 0, sd = unit_fe_sd)
  baseline_logit <- stats::qlogis(base_prob)

  # Covariates
  if (add_covariates) {
    x1 <- stats::rnorm(n)
    x2 <- stats::rbinom(n, 1, 0.5)
    # Make x1 predictive of both treatment and outcome
    alpha_i <- alpha_i + 0.3 * x1
    group_prob_adj <- group_vec  # keep groups fixed for simplicity
  }

  # Build panel
  rows <- vector("list", n * nperiods)
  k    <- 0L

  for (i in seq_len(n)) {
    g_i <- group_vec[i]

    for (t in seq_len(nperiods)) {
      k <- k + 1L
      treated <- (g_i > 0 && t >= g_i)
      att_g   <- if (g_i > 0) cohort_att_map[as.character(g_i)] else 0
      att_g   <- if (is.na(att_g)) 0 else att_g

      # Time FE (gradual trend)
      lambda_t <- 0.05 * (t - 1)

      # Latent utility
      epsilon <- stats::rlogis(1)
      eta     <- baseline_logit + alpha_i[i] + lambda_t +
        att_g * treated + epsilon

      y_it <- as.integer(eta > 0)

      row_data <- list(
        id       = i,
        period   = t,
        y        = y_it,
        g        = g_i,
        D        = as.integer(treated),
        alpha_i  = round(alpha_i[i], 4)
      )

      if (add_covariates) {
        row_data$x1 <- round(x1[i], 4)
        row_data$x2 <- x2[i]
      }

      rows[[k]] <- row_data
    }
  }

  dat <- do.call(rbind, lapply(rows, as.data.frame))
  rownames(dat) <- NULL

  attr(dat, "true_att")        <- true_att
  attr(dat, "cohort_periods")  <- cohort_periods
  attr(dat, "dgp")             <- "binary_logit"

  dat
}


#' @title Simulate Count Panel Data with Staggered Treatment
#'
#' @description
#' Generates simulated panel data with a count outcome (Poisson-distributed)
#' and staggered treatment adoption. Treatment effect is multiplicative
#' (rate ratio) on the count scale.
#'
#' @param n Integer. Number of units. Default 500.
#' @param nperiods Integer. Number of time periods. Default 6.
#' @param prop_treated Numeric. Proportion of units ever treated. Default 0.5.
#' @param n_cohorts Integer. Number of treatment cohorts. Default 3.
#' @param true_rr Numeric or vector. True rate ratio for each cohort.
#'   Default 1.5 (50 percent increase in count).
#' @param base_rate Numeric. Baseline Poisson rate. Default 5.
#' @param overdispersion Logical. Add overdispersion (negative binomial).
#'   Default FALSE.
#' @param seed Integer. Random seed.
#'
#' @return Long-format data frame with columns: id, period, y, g, D, x1.
#'
#' @examples
#' dat <- sim_count_panel(n = 400, nperiods = 6, true_rr = 1.8)
#' summary(dat$y)
#'
#' @export
sim_count_panel <- function(
    n              = 500L,
    nperiods       = 6L,
    prop_treated   = 0.5,
    n_cohorts      = 3L,
    true_rr        = 1.5,
    base_rate      = 5,
    overdispersion = FALSE,
    seed           = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  if (length(true_rr) == 1) true_rr <- rep(true_rr, n_cohorts)

  n_treat      <- round(n * prop_treated)
  n_never      <- n - n_treat
  n_per_cohort <- floor(n_treat / n_cohorts)
  cohort_periods <- round(seq(2, nperiods - 1, length.out = n_cohorts))

  group_vec <- c(rep(0, n_never),
                 unlist(lapply(seq_len(n_cohorts), function(k)
                   rep(cohort_periods[k], n_per_cohort))))
  while (length(group_vec) < n) group_vec <- c(group_vec, 0L)
  group_vec <- group_vec[seq_len(n)]

  rr_map  <- stats::setNames(true_rr, cohort_periods)
  alpha_i <- exp(stats::rnorm(n, mean = 0, sd = 0.3))
  x1      <- stats::rnorm(n)

  rows <- vector("list", n * nperiods)
  k    <- 0L

  for (i in seq_len(n)) {
    g_i <- group_vec[i]
    for (t in seq_len(nperiods)) {
      k <- k + 1L
      treated <- (g_i > 0 && t >= g_i)
      rr_g    <- if (g_i > 0) rr_map[as.character(g_i)] else 1
      rr_g    <- if (is.na(rr_g)) 1 else rr_g

      lambda_t <- 1 + 0.1 * (t - 1)
      mu_it    <- base_rate * alpha_i[i] * lambda_t * (rr_g ^ treated)

      y_it <- if (overdispersion) {
        MASS::rnegbin(1, mu = mu_it, theta = 5)
      } else {
        stats::rpois(1, lambda = mu_it)
      }

      rows[[k]] <- data.frame(id = i, period = t, y = y_it,
                              g = g_i, D = as.integer(treated),
                              x1 = round(x1[i], 4))
    }
  }

  dat <- do.call(rbind, rows)
  rownames(dat) <- NULL
  attr(dat, "true_rr")        <- true_rr
  attr(dat, "cohort_periods") <- cohort_periods
  attr(dat, "dgp")            <- "count_poisson"
  dat
}
