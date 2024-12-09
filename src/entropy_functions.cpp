// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// ----------------------------------------------------------------------------
// Function: entropySingle
// Description:
//   Calculates the entropy of a probability vector using the formula:
//     Entropy = -sum(x_i * log2(x_i)) for all x_i > 0
// Parameters:
//   p - Eigen::VectorXd representing a probability distribution (non-negative and sums to 1)
// Returns:
//   double representing the entropy of the distribution
// ----------------------------------------------------------------------------

// [[Rcpp::export]]
double entropySingle(const Eigen::VectorXd& p) {
  // Input Validation: Ensure that all probabilities are non-negative
  if ( (p.array() < 0.0).any() ) {
    Rcpp::stop("Input vector 'p' contains negative probabilities.");
  }

  // Optional Input Validation: Check if the vector sums to 1 within a tolerance
  double sum_p = p.sum();
  if ( std::abs(sum_p - 1.0) > 1e-6 ) {
    Rcpp::warning("Input vector 'p' does not sum to 1 (sum = %.6f). Proceeding with entropy calculation.", sum_p);
  }

  // Compute entropy terms using Eigen's vectorization with a lambda function
  // For each element x_i in p:
  //   If x_i > 0, compute -x_i * log2(x_i)
  //   Else, assign 0 (since 0 * log2(0) is defined as 0)
  Eigen::VectorXd entropyTerms = p.unaryExpr([](double x) -> double {
    return (x > 0.0) ? (-x * std::log2(x)) : 0.0;
  });

  // Sum the entropy terms to get the total entropy
  double entropy = entropyTerms.sum();

  return entropy;
}
