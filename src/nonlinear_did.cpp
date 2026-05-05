// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

//' @title Fast multiplier bootstrap weights
//' @description Generates exponential multiplier weights for the
//'   wild/multiplier bootstrap used in nonlinear DiD inference.
//' @param n Integer. Number of units.
//' @param nboot Integer. Number of bootstrap iterations.
//' @return A matrix of dimension n x nboot with positive weights
//'   that sum to n in each column.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix multiplier_weights_cpp(int n, int nboot) {
  NumericMatrix weights(n, nboot);
  for (int b = 0; b < nboot; b++) {
    NumericVector w = rexp(n, 1.0);
    double mean_w = mean(w);
    for (int i = 0; i < n; i++) {
      weights(i, b) = w[i] / mean_w;
    }
  }
  return weights;
}


//' @title Weighted mean (fast)
//' @description Computes weighted mean efficiently for bootstrap.
//' @param x Numeric vector.
//' @param w Numeric weights vector.
//' @return Weighted mean scalar.
//' @keywords internal
// [[Rcpp::export]]
double weighted_mean_cpp(NumericVector x, NumericVector w) {
  int n = x.size();
  double num = 0.0, den = 0.0;
  for (int i = 0; i < n; i++) {
    if (!NumericVector::is_na(x[i]) && !NumericVector::is_na(w[i])) {
      num += w[i] * x[i];
      den += w[i];
    }
  }
  if (den == 0.0) return NA_REAL;
  return num / den;
}


//' @title Logit-scale DR score computation
//' @description Computes the doubly-robust score function for the
//'   logit-scale DiD estimator. Used internally for fast SE computation.
//' @param delta_logit Numeric vector of logit(Y1) - logit(Y0).
//' @param D Integer vector of treatment indicators.
//' @param pscore Numeric vector of propensity scores.
//' @param mu_hat Numeric vector of outcome model predictions.
//' @return Numeric vector of DR scores (influence functions).
//' @keywords internal
// [[Rcpp::export]]
NumericVector dr_score_logit_cpp(NumericVector delta_logit,
                                  IntegerVector D,
                                  NumericVector pscore,
                                  NumericVector mu_hat) {
  int n = delta_logit.size();
  NumericVector score(n);
  double p_D = mean(as<NumericVector>(D));
  if (p_D <= 0.0) return score;

  for (int i = 0; i < n; i++) {
    double w_t = D[i] / p_D;
    double w_c = (1.0 - D[i]) * pscore[i] / (1.0 - pscore[i]) / p_D;
    score[i] = w_t * (delta_logit[i] - mu_hat[i]) -
               w_c * (delta_logit[i] - mu_hat[i]);
  }
  return score;
}


//' @title Compute ATT from DR scores
//' @description Returns ATT as mean of DR influence function scores.
//' @param score Numeric vector from dr_score_logit_cpp.
//' @return Scalar ATT estimate.
//' @keywords internal
// [[Rcpp::export]]
double att_from_score_cpp(NumericVector score) {
  return mean(score);
}
