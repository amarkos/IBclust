// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Jensen–Shannon divergence between two row‐vectors
double js_div(const arma::rowvec &p, const arma::rowvec &q) {
  arma::rowvec m = 0.5 * (p + q);
  arma::rowvec pm = arma::log2(p / m + 1e-16);
  arma::rowvec qm = arma::log2(q / m + 1e-16);
  return 0.5 * (arma::dot(p, pm) + arma::dot(q, qm));
}

// [[Rcpp::export]]
arma::mat make_IB_distmat(const arma::mat &p_xy) {
  int n = p_xy.n_rows;
  arma::vec pz = arma::sum(p_xy, 1);             // p(z)
  arma::mat p_ygz = p_xy.each_col() / pz;        // p(y|z)
  arma::mat D(n,n, arma::fill::zeros);

  for(int i = 0; i < n; ++i) {
    for(int j = i+1; j < n; ++j) {
      double cost = (pz[i] + pz[j]) * js_div(p_ygz.row(i), p_ygz.row(j));
      D(i,j) = D(j,i) = cost;
    }
  }
  return D;
}

// [[Rcpp::export]]
double js_divergence(const NumericVector &p, const NumericVector &q) {
  // convenience wrapper to call from R
  arma::rowvec pv = as<arma::rowvec>(p);
  arma::rowvec qv = as<arma::rowvec>(q);
  return js_div(pv, qv);
}

// [[Rcpp::export]]
double mutual_information(const arma::mat& p_xy, double base = 2.0) {
  // Marginals
  arma::vec p_x = arma::sum(p_xy, 1);            // row sums
  arma::vec p_y = arma::sum(p_xy, 0).t();        // column sums

  double mi = 0.0;
  double log_base = std::log(base);

  arma::uword R = p_xy.n_rows;
  arma::uword C = p_xy.n_cols;
  for (arma::uword i = 0; i < R; ++i) {
    for (arma::uword j = 0; j < C; ++j) {
      double p = p_xy(i, j);
      if (p > 0) {
        // log(p / (p_x * p_y)) in given base
        mi += p * (std::log(p / (p_x(i) * p_y(j))) / log_base);
      }
    }
  }

  return mi;
}
