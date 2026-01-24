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

// -----------------------------------------------------------------------------
// Compute JS divergence between two distributions from Nystrom factors
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double js_divergence_nystrom(const arma::vec& p_z_i, const arma::vec& p_z_j) {
  arma::vec m = 0.5 * (p_z_i + p_z_j);
  
  double kl_p_m = 0.0;
  double kl_q_m = 0.0;
  
  for (unsigned int k = 0; k < m.n_elem; k++) {
    if (p_z_i(k) > 1e-300 && m(k) > 1e-300) {
      kl_p_m += p_z_i(k) * std::log2(p_z_i(k) / m(k));
    }
    if (p_z_j(k) > 1e-300 && m(k) > 1e-300) {
      kl_q_m += p_z_j(k) * std::log2(p_z_j(k) / m(k));
    }
  }
  
  return 0.5 * kl_p_m + 0.5 * kl_q_m;
}

// -----------------------------------------------------------------------------
// Compute IB distance matrix using Nystrom approximation
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat make_IB_distmat_nystrom_cpp(const arma::mat& B,
                                       const arma::vec& col_sums) {
  int n = B.n_rows;
  arma::mat D(n, n, arma::fill::zeros);
  
  // Compute p(y|x) for each observation
  
  // First compute the unnormalised kernel weights for each obs
  // Precompute all conditional distributions p(y|x)
  std::vector<arma::vec> p_y_x(n);
  arma::vec row_probs(n); 
  
  for (int i = 0; i < n; i++) {
    arma::vec unnorm = B * B.row(i).t();  // n x 1
    
    // Normalise by col_sums
    if (col_sums(i) > 1e-10) {
      unnorm = unnorm / col_sums(i);
    }
    
    // Normalise
    double total = arma::accu(unnorm);
    if (total > 1e-10) {
      p_y_x[i] = unnorm / total;
      row_probs(i) = total;
    } else {
      p_y_x[i] = arma::vec(n, arma::fill::zeros);
      p_y_x[i](i) = 1.0; 
      row_probs(i) = 1.0 / n;
    }
  }
  
  row_probs = row_probs / arma::accu(row_probs);
  
  // pairwise JS divergences
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      double js = js_divergence_nystrom(p_y_x[i], p_y_x[j]);
      double weight = row_probs(i) + row_probs(j);
      D(i, j) = weight * js;
      D(j, i) = D(i, j);
    }
  }
  
  return D;
}

// -----------------------------------------------------------------------------
// Compute mutual information from Nystrom factors
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
double mutual_information_nystrom_cpp(const arma::mat& B,
                                       const arma::vec& col_sums) {
  int n = B.n_rows;
  
  double px = 1.0 / n;
  
  // Compute p(y) = sum_x p(x,y) = sum_x p(x) * p(y|x)
  arma::vec p_y(n, arma::fill::zeros);
  std::vector<arma::vec> p_y_x(n);
  
  for (int i = 0; i < n; i++) {
    arma::vec unnorm = B * B.row(i).t();
    if (col_sums(i) > 1e-10) {
      unnorm = unnorm / col_sums(i);
    }
    double total = arma::accu(unnorm);
    if (total > 1e-10) {
      p_y_x[i] = unnorm / total;
    } else {
      p_y_x[i] = arma::vec(n, arma::fill::zeros);
      p_y_x[i](i) = 1.0;
    }
    p_y += px * p_y_x[i];
  }
  
  // Compute MI = sum_x,y p(x,y) * log2(p(y|x) / p(y))
  double mi = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double p_xy = px * p_y_x[i](j);
      if (p_xy > 1e-300 && p_y(j) > 1e-300 && p_y_x[i](j) > 1e-300) {
        mi += p_xy * std::log2(p_y_x[i](j) / p_y(j));
      }
    }
  }
  
  return mi;
}
