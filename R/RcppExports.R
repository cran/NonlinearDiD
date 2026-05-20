# Pure-R implementations replacing the former compiled Rcpp helpers.
# These are functionally identical; the C++ versions were speed
# optimisations only. Removing the compiled dependency simplifies
# installation and CRAN submission.

#' @keywords internal
multiplier_weights_cpp <- function(n, nboot) {
  # n x nboot matrix of Exp(1) weights, each column normalised to sum to n
  m <- matrix(stats::rexp(n * nboot), nrow = n, ncol = nboot)
  sweep(m, 2, colMeans(m), "/")
}

#' @keywords internal
weighted_mean_cpp <- function(x, w) {
  sum(x * w) / sum(w)
}

#' @keywords internal
dr_score_logit_cpp <- function(delta_logit, D, pscore, mu_hat) {
  p_D <- mean(D)
  w_t <- D / p_D
  w_c <- (1 - D) * pscore / (1 - pscore) / p_D
  w_t * (delta_logit - mu_hat) - w_c * (delta_logit - mu_hat)
}

#' @keywords internal
att_from_score_cpp <- function(score) {
  mean(score)
}

