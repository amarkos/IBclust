#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// -----------------------------------------------------------------------------
// Function to calculate KL divergence for two vectors
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
double klSingle(const arma::vec& p, const arma::vec& q) {
  // Ensure p and q have the same length
  if (p.n_elem != q.n_elem) {
    Rcpp::stop("Vectors 'p' and 'q' must have the same length.");
  }

  // Identify indices where both p and q are greater than zero
  arma::uvec valid = arma::find((p > 0) && (q > 0));

  // If no valid elements, return zero and issue a warning
  if(valid.n_elem == 0){
    Rcpp::warning("No valid elements to compute KL divergence. Returning 0.");
    return 0.0;
  }

  // Calculate KL divergence only for valid elements
  arma::vec kl_terms = p.elem(valid) % arma::log(p.elem(valid) / q.elem(valid));

  return arma::accu(kl_terms);
}

// -----------------------------------------------------------------------------
// Function to compute the natural logarithm of a vector with validation
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
arma::vec vlog(const arma::vec& x) {
  // Compute the logarithm
  arma::vec log_x = arma::log(x);

  // Identify non-positive elements
  arma::uvec non_positive = arma::find(x <= 0);
  if(non_positive.n_elem > 0){
    Rcpp::warning("Non-positive elements detected in input to 'vlog'. Assigning -Inf to these positions.");
    log_x.elem(non_positive).fill(-arma::datum::inf);
  }

  return log_x;
}

// -----------------------------------------------------------------------------
// Function to perform the qt_x_step in C++
// Assigns each data point to the cluster with the highest adjusted KL divergence
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat qt_x_step_cpp(int n_rows, int T, double beta,
                        const arma::mat& py_x, const arma::mat& qy_t, const arma::vec& qt) {

  // Initialize qt_x matrix with zeros (Clusters x Data Points)
  arma::mat qt_x(T, n_rows, arma::fill::zeros);

  // Precompute log(qt) for efficiency
  arma::vec log_qt = vlog(qt);

  // Iterate over each data point
  for (int x = 0; x < n_rows; ++x) {
    arma::vec kl_divs(T);

    // Compute KL divergence between data point x and each cluster
    for (int t = 0; t < T; ++t) {
      kl_divs(t) = klSingle(py_x.col(x), qy_t.col(t));
    }

    // Compute the adjusted divergence
    arma::vec l = log_qt - beta * kl_divs;

    // Identify the cluster with the maximum adjusted divergence
    arma::uword t_max = l.index_max(); // FIXED: Use index_max() instead of max()
    
    // Assign the data point to the identified cluster
    qt_x(t_max, x) = 1;
  }

  return qt_x;
}

// -----------------------------------------------------------------------------
// Main function to adjust beta and update cluster assignments
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
List qt_x_step_beta_cpp(int n_rows, int T,
                        const arma::mat& py_x, const arma::mat& qy_t,
                        const arma::vec& qt, arma::mat qt_x) {

  // Find the cluster with the minimum value in qt
  arma::uword min_t = qt.index_min();

  // Identify all data points assigned to the cluster 'min_t'
  arma::uvec x_min_t = arma::find(qt_x.row(min_t) == 1);

  // Handle the case with a single cluster
  if (T == 1) {
    qt_x.fill(1);
    return List::create(Named("qt_x") = qt_x);
  } else {
    // Initialize a vector to store minimum beta values for each relevant data point
    arma::vec beta_min(x_min_t.n_elem, arma::fill::zeros);

    // Precompute log(qt) for efficiency
    arma::vec log_qt = vlog(qt);

    // Iterate over each relevant data point
    for (arma::uword i = 0; i < x_min_t.n_elem; ++i) {
      int x = x_min_t[i];
      arma::vec kl_divs(T);

      // Compute KL divergence between data point x and each cluster
      for (int t = 0; t < T; ++t) {
        kl_divs(t) = klSingle(py_x.col(x), qy_t.col(t));
      }

      // Compute the ratio for beta adjustment
      arma::vec numerator = log_qt[min_t] - log_qt;
      arma::vec denominator = kl_divs[min_t] - kl_divs;

      // Prevent division by zero by replacing zero denominators with a small value
      arma::uvec zero_denominator = arma::find(denominator == 0);
      if(zero_denominator.n_elem > 0){
        denominator.elem(zero_denominator).fill(1e-10);
      }

      arma::vec ratio = numerator / denominator;

      // Assign the maximum ratio to beta_min
      beta_min(i) = arma::max(ratio);
    }

    // Determine the new beta value
    double beta = std::max(arma::max(beta_min) + 1e-3, 1.0);

    // Update qt_x with the new beta
    qt_x = qt_x_step_cpp(n_rows, T, beta, py_x, qy_t, qt);

    return List::create(Named("qt_x") = qt_x, Named("beta") = beta);
  }
}

// -----------------------------------------------------------------------------
// Function to update qy_t matrix based on current assignments and probabilities
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat qy_t_step_cpp(const arma::mat& py_x, const arma::mat& qt_x,
                        const arma::vec& qt, const arma::vec& px) {
  // Check for zeros in qt to prevent division by zero
  arma::uvec zero_qt = arma::find(qt == 0);
  arma::vec safe_qt = qt;
  if(zero_qt.n_elem > 0){
    Rcpp::warning("Zero elements detected in 'qt'. Assigning a small value to prevent division by zero.");
    safe_qt.elem(zero_qt).fill(1e-10);
  }

  // Compute the inverse of qt
  arma::vec inv_qt = 1.0 / safe_qt;

  // Perform the outer product (1/qt) * px.t(), resulting in a (T x n) matrix
  arma::mat qt_px = inv_qt * px.t();

  // Element-wise multiply qt_x by qt_px
  arma::mat scaled_qt_x = qt_x % qt_px;

  // Compute the final qy_t matrix using matrix multiplication
  arma::mat qy_t = py_x * scaled_qt_x.t();

  return qy_t;
}

// -----------------------------------------------------------------------------
// Function to perform the qt_x_step for the IB in C++
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat qt_x_step_gib_cpp(int n_rows, int T, double beta, double alpha,
                           const arma::mat& py_x, const arma::mat& qy_t, const arma::vec& qt) {

  // Initialize qt_x matrix with zeros (Clusters x Data Points)
  arma::mat qt_x(T, n_rows, arma::fill::zeros);

  // Precompute log(qt) for efficiency
  arma::vec log_qt = vlog(qt);

  // Iterate over each data point
  for (int x = 0; x < n_rows; ++x) {
    arma::vec kl_divs(T);

    // Compute KL divergence between data point x and each cluster
    for (int t = 0; t < T; ++t) {
      kl_divs(t) = klSingle(py_x.col(x), qy_t.col(t));
    }

    // Compute the adjusted divergence
    arma::vec l = arma::exp( (log_qt - beta * kl_divs) / alpha);

    // Normalise each column to sum to 1
    double S = arma::accu(l);
    if (S > 0.0) {
      qt_x.col(x) = l / S;
    } else {
      // fallback to uniform if underflow/zero
      qt_x.col(x).fill(1.0 / T);
    }
  }

  return qt_x;
}
