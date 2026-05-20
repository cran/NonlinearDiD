#' NonlinearDiD: Staggered DiD with Nonlinear Outcomes
#'
#' @description
#' \strong{NonlinearDiD} supports staggered difference-in-differences designs
#' with nonlinear outcomes for both \emph{panel} and \emph{repeated
#' cross-section} data.
#'
#' For panel data, units are followed over time and \code{idname} identifies
#' repeated observations. For repeated cross-section data, observations are
#' independent within each time period; \code{idname} is optional and may
#' identify survey records or households, but the estimator does not require
#' the same units to appear across periods.
#'
#' The package extends the Callaway and Sant'Anna (2021) framework to
#' nonlinear outcome models, including binary (logit/probit), count
#' (Poisson/NegBin), and odds-ratio estimands.
#'
#' \strong{The Core Problem}
#'
#' The canonical CS2021 framework assumes parallel trends on the mean scale
#' of a continuous outcome. For binary and count outcomes, this assumption
#' is not scale-invariant: parallel trends in P(Y=1) does NOT imply parallel
#' trends in log-odds, pre-trend tests depend on which scale is used, and
#' treatment effect estimates conflate true effects with Jensen's inequality.
#'
#' \strong{Main Functions}
#'
#' \itemize{
#'   \item \code{nonlinear_attgt()} -- Estimate ATT(g,t) under nonlinear
#'     outcome models; supports panel and repeated cross-sections, with
#'     optional sampling weights (\code{weightsname}) and clustered
#'     inference (\code{cluster_var}).
#'   \item \code{nonlinear_aggte()} -- Aggregate: event-study, group,
#'     calendar, overall.
#'   \item \code{nonlinear_pretest()} -- Pre-treatment parallel trends test.
#'   \item \code{binary_did_logit()} -- 2x2 DiD with logit outcome.
#'   \item \code{binary_did_probit()} -- 2x2 DiD with probit outcome.
#'   \item \code{binary_did_dr()} -- Doubly-robust binary DiD.
#'   \item \code{count_did_poisson()} -- Poisson QMLE DiD for count outcomes.
#'   \item \code{odds_ratio_did()} -- Odds-ratio DiD (scale-free).
#'   \item \code{nonlinear_bounds()} -- Nonparametric Manski / PT bounds.
#'   \item \code{sim_binary_panel()} -- Simulate binary staggered panel data.
#'   \item \code{sim_count_panel()} -- Simulate count staggered panel data.
#'   \item \code{sim_binary_rcs()} -- Simulate binary repeated cross-section
#'     data.
#' }
#'
#' \strong{Quick Start: Panel}
#'
#' \preformatted{
#' library(NonlinearDiD)
#' dat <- sim_binary_panel(n = 500, nperiods = 8, seed = 42)
#' res <- nonlinear_attgt(dat, yname = "y", tname = "period",
#'                         idname = "id", gname = "g",
#'                         outcome_model = "logit")
#' agg <- nonlinear_aggte(res, type = "dynamic")
#' plot(agg)
#' nonlinear_pretest(res)
#' }
#'
#' \strong{Quick Start: Repeated Cross-Section}
#'
#' \preformatted{
#' library(NonlinearDiD)
#' rcs <- sim_binary_rcs(n_per_period = 500, nperiods = 8, seed = 7)
#' res <- nonlinear_attgt(rcs, yname = "y", tname = "period",
#'                         gname = "g", outcome_model = "logit",
#'                         data_type = "repeated_cross_section",
#'                         estimand = "ape",
#'                         control_group = "notyetreated")
#' plot(nonlinear_aggte(res, type = "dynamic"))
#' }
#'
#' \strong{Survey-Weighted Repeated Cross-Section Example}
#'
#' \preformatted{
#' # Example: CPS-FSS-style data with survey weights and state clustering
#' # res <- nonlinear_attgt(
#' #   data          = my_survey_data,
#' #   yname         = "food_insecure",
#' #   tname         = "year",
#' #   gname         = "policy_end_year",
#' #   idname        = "household_id",
#' #   data_type     = "repeated_cross_section",
#' #   outcome_model = "logit",
#' #   estimand      = "ape",
#' #   weightsname   = "survey_weight",
#' #   cluster_var   = "state",
#' #   control_group = "notyetreated"
#' # )
#' }
#'
#' @references
#' Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with
#' multiple time periods. \emph{Journal of Econometrics}, 225(2), 200-230.
#'
#' Roth, J., & Sant'Anna, P. H. C. (2023). When is parallel trends sensitive
#' to functional form? \emph{Econometrica}, 91(2), 737-747.
#'
#' Wooldridge, J. M. (2023). Simple approaches to nonlinear
#' difference-in-differences with panel data. \emph{The Econometrics Journal}, 26(3).
#'
#' Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust
#' difference-in-differences estimators. \emph{Journal of Econometrics},
#' 219(1), 101-122.
#'
#' @docType package
#' @name NonlinearDiD-package
#' @aliases NonlinearDiD
"_PACKAGE"
