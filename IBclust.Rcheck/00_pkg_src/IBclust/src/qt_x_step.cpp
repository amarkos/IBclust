#include <RcppArmadillo.h>
#include <cmath>
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


// [[Rcpp::export]]
List qt_x_step_cpp(int n_rows, int T, double beta,
                   const arma::mat& py_x, const arma::mat& qy_t, const arma::vec& qt) {

  // Initialize qt_x matrix with zeros (Clusters x Data Points)
  arma::mat qt_x(T, n_rows, arma::fill::zeros);

  // Precompute log(qt) for efficiency
  arma::vec log_qt = vlog(qt);

  // Store the loss matrix for potential cluster rescue
  arma::mat L(T, n_rows);

  // Iterate over each data point
  for (int x = 0; x < n_rows; ++x) {
    arma::vec kl_divs(T);

    // Compute KL divergence between data point x and each cluster
    for (int t = 0; t < T; ++t) {
      kl_divs(t) = klSingle(py_x.col(x), qy_t.col(t));
    }

    // Compute the adjusted divergence
    arma::vec l = log_qt - beta * kl_divs;
    L.col(x) = l;  // Store for later

    // Identify the cluster with the maximum adjusted divergence
    arma::uword t_max = l.index_max();
    
    // Assign the data point to the identified cluster
    qt_x(t_max, x) = 1;
  }

  // Check for empty clusters and rescue
  arma::vec cluster_sizes = arma::sum(qt_x, 1);
  arma::uvec empty_clusters = arma::find(cluster_sizes == 0);
  
  bool rescue_occurred = false;  // Track if rescue happened
  
  if (empty_clusters.n_elem > 0) {
    rescue_occurred = true;  // Mark that rescue was needed
    
    // For each empty cluster
    for (arma::uword i = 0; i < empty_clusters.n_elem; ++i) {
      arma::uword empty_c = empty_clusters(i);
      
      // Find observations in clusters with size > 1
      arma::uvec current_assignments(n_rows);
      for (int x = 0; x < n_rows; ++x) {
        current_assignments(x) = arma::index_max(qt_x.col(x));
      }
      
      // Build list of observations that can be borrowed
      std::vector<arma::uword> available_obs;
      for (int x = 0; x < n_rows; ++x) {
        arma::uword curr_cluster = current_assignments(x);
        if (cluster_sizes(curr_cluster) > 1) {
          available_obs.push_back(x);
        }
      }
      
      // If no observations can be safely borrowed, steal from any cluster
      if (available_obs.empty()) {
        for (int x = 0; x < n_rows; ++x) {
          available_obs.push_back(x);
        }
      }
      
      // Find the observation with the best loss for the empty cluster
      double best_loss = -arma::datum::inf;
      arma::uword best_obs = available_obs[0];
      
      for (arma::uword x : available_obs) {
        if (L(empty_c, x) > best_loss) {
          best_loss = L(empty_c, x);
          best_obs = x;
        }
      }
      
      // Borrow this observation
      arma::uword old_cluster = current_assignments(best_obs);
      qt_x(old_cluster, best_obs) = 0;
      qt_x(empty_c, best_obs) = 1;
      
      // Update cluster sizes
      cluster_sizes(old_cluster) -= 1;
      cluster_sizes(empty_c) += 1;
    }
  }

  // Return both the matrix and the rescue flag
  return List::create(
    Named("qt_x") = qt_x,
    Named("rescue_occurred") = rescue_occurred
  );
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
    return List::create(
      Named("qt_x") = qt_x,
      Named("beta") = 1.0,
      Named("rescue_occurred") = false
    );
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

    // Update qt_x with the new beta - NOW RETURNS A LIST
    List qt_x_result = qt_x_step_cpp(n_rows, T, beta, py_x, qy_t, qt);
    
    // Extract the matrix and rescue flag from the returned List
    arma::mat qt_x_new = qt_x_result["qt_x"];
    bool rescue_occurred = qt_x_result["rescue_occurred"];

    return List::create(
      Named("qt_x") = qt_x_new, 
      Named("beta") = beta,
      Named("rescue_occurred") = rescue_occurred
    );
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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat qy_t_step_nystrom_cpp(const arma::mat& B,
                                 const arma::vec& col_sums,
                                 const arma::mat& qt_x,
                                 const arma::vec& qt,
                                 const arma::vec& px) {
  
  int n = B.n_rows;
  
  // Check for zeros in qt to prevent division by zero
  arma::vec safe_qt = qt;
  arma::uvec zero_qt = arma::find(qt == 0);
  if(zero_qt.n_elem > 0) {
    Rcpp::warning("Zero elements detected in 'qt'. Assigning small value.");
    safe_qt.elem(zero_qt).fill(1e-10);
  }
  
  // Compute inv_qt and qt_px
  arma::vec inv_qt = 1.0 / safe_qt;
  arma::mat qt_px = inv_qt * px.t(); 
  
  // Element-wise multiplication of qt_x by qt_px
  arma::mat scaled_qt_x = qt_x % qt_px;
 
  
  // Compute D * scaled_qt_x.t() 
  arma::mat scaled_qt_x_t = scaled_qt_x.t(); 
  
  for(int i = 0; i < n; i++) {
    if(col_sums(i) > 1e-10) { 
      scaled_qt_x_t.row(i) /= col_sums(i);
    } else {
      scaled_qt_x_t.row(i).zeros();
    }
  }
  
  // Compute B * (B^T * scaled_qt_x_t)
  arma::mat temp = B.t() * scaled_qt_x_t;  
  arma::mat qy_t = B * temp; 
  
  return qy_t;
}

// [[Rcpp::depends(RcppArmadillo)]]
// -----------------------------------------------------------------------------
// HELPER: Compute Cross-Entropy Matrix efficiently
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat compute_cross_entropy_nystrom(const arma::mat& B,
                                        const arma::vec& col_sums,
                                        const arma::mat& qy_t) {
  
  // Log of centroids (safe)
  arma::mat log_qy_t = arma::log(qy_t);
  
  // Use small value -1e100 instead of -Inf to avoid NaN propagation
  log_qy_t.elem(arma::find(qy_t <= 1e-300)).fill(-1e100); 
  
  // Efficient Chain Multiplication
  arma::mat temp = B.t() * log_qy_t;
  arma::mat M_raw = B * temp;
  
  // Row normalisation
  arma::mat M = M_raw;
  for(unsigned int i = 0; i < M.n_rows; ++i) {
    if(col_sums(i) > 1e-10) {
      M.row(i) /= col_sums(i);
    } else {
      M.row(i).zeros();
    }
  }
  
  return M;
}

// -----------------------------------------------------------------------------
// FAST REPLACEMENT: qt_x_step_nystrom_cpp
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List qt_x_step_nystrom_cpp(int n_rows, int T, double beta,
                           const arma::mat& B,
                           const arma::vec& col_sums,
                           const arma::mat& qy_t,
                           const arma::vec& qt) {

  arma::mat qt_x(T, n_rows, arma::fill::zeros);
  arma::vec log_qt = arma::log(qt);
  log_qt.elem(arma::find(qt <= 1e-300)).fill(-1e100);

  // Precompute Cross-Entropy Matrix M
  arma::mat M = compute_cross_entropy_nystrom(B, col_sums, qy_t);
  
  // Iterative procedure
  // Set score = log q(t) + beta * M(x, t)
  for (int x = 0; x < n_rows; ++x) {
    arma::vec scores = log_qt + beta * M.row(x).t();
    arma::uword t_max = scores.index_max();
    qt_x(t_max, x) = 1;
  }

  // Rescue Empty Clusters
  arma::vec cluster_sizes = arma::sum(qt_x, 1);
  arma::uvec empty_clusters = arma::find(cluster_sizes == 0);
  bool rescue_occurred = false;
  
  if (empty_clusters.n_elem > 0) {
    rescue_occurred = true;
    
    // Current assignments
    arma::uvec current_assignments(n_rows);
    for (int x = 0; x < n_rows; ++x) {
      current_assignments(x) = arma::index_max(qt_x.col(x));
    }
    
    for (arma::uword i = 0; i < empty_clusters.n_elem; ++i) {
      arma::uword empty_c = empty_clusters(i);
      
      // Candidate pool from clusters with size > 1
      std::vector<arma::uword> available_obs;
      for (int x = 0; x < n_rows; ++x) {
        if (cluster_sizes(current_assignments(x)) > 1) {
          available_obs.push_back(x);
        }
      }
      
      // Fallback: if all clusters size 1, allow stealing
      if (available_obs.empty()) {
        for(int x=0; x<n_rows; ++x) available_obs.push_back(x);
      }
      
      // Pick observation with max cross entropy
      double best_val = -arma::datum::inf;
      arma::uword best_obs = available_obs[0];
      
      for (arma::uword x : available_obs) {
        // Safe access
        double val = M(x, empty_c);
        if (val > best_val) {
          best_val = val;
          best_obs = x;
        }
      }
      
      arma::uword old_cluster = current_assignments(best_obs);
      qt_x(old_cluster, best_obs) = 0;
      qt_x(empty_c, best_obs) = 1;
      
      // Update sizes locally to prevent stealing same point twice
      cluster_sizes(old_cluster)--;
      cluster_sizes(empty_c)++;
      current_assignments(best_obs) = empty_c;
    }
  }

  return List::create(
    Named("qt_x") = qt_x,
    Named("rescue_occurred") = rescue_occurred
  );
}

// -----------------------------------------------------------------------------
// FAST REPLACEMENT: qt_x_step_beta_nystrom_cpp
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List qt_x_step_beta_nystrom_cpp(int n_rows, int T,
                                const arma::mat& B,
                                const arma::vec& col_sums,
                                const arma::mat& qy_t,
                                const arma::vec& qt,
                                const arma::mat& qt_x) {

  arma::uword min_t = qt.index_min();
  arma::uvec x_min_t = arma::find(qt_x.row(min_t) == 1);

  // Safety check for single point cluster
  if (T == 1) {
    arma::mat qt_x_out = qt_x; 
    qt_x_out.fill(1);
    return List::create(Named("qt_x") = qt_x_out, Named("beta") = 1.0, Named("rescue_occurred") = false);
  } 
  
  // Sefety check for empty cluster case
  // Return beta = 10.0 and let assignment step rescue cluster
  if (x_min_t.n_elem == 0) {
    double default_beta = 10.0;
    List qt_x_result = qt_x_step_nystrom_cpp(n_rows, T, default_beta, B, col_sums, qy_t, qt);
    
    // Pass through 
    return List::create(
      Named("qt_x") = qt_x_result["qt_x"],
      Named("beta") = default_beta,
      Named("rescue_occurred") = true 
    );
  }

  // Precompute Cross-Entropy Matrix M
  arma::mat M = compute_cross_entropy_nystrom(B, col_sums, qy_t);
  
  arma::vec beta_min(x_min_t.n_elem, arma::fill::zeros);
  arma::vec log_qt = arma::log(qt);
  log_qt.elem(arma::find(qt <= 1e-300)).fill(-1e100);

  // Beta constraints
  for (arma::uword i = 0; i < x_min_t.n_elem; ++i) {
    int x = x_min_t[i];
    
    arma::vec denom = M.row(x).t() - M(x, min_t);
    arma::vec num = log_qt(min_t) - log_qt;
    
    arma::uvec zero_denom = arma::find(arma::abs(denom) < 1e-9);
    if (zero_denom.n_elem > 0) denom.elem(zero_denom).fill(1e-10);
    
    arma::vec ratios = num / denom;
    
    // Safety check: filter out NaNs or Infs if any
    double max_r = -arma::datum::inf;
    for(unsigned int k=0; k<ratios.n_elem; k++) {
        if(std::isfinite(ratios(k)) && ratios(k) > max_r) {
            max_r = ratios(k);
        }
    }
    // No valid ratio found => default to 0
    if(max_r == -arma::datum::inf) max_r = 0.0;
    
    beta_min(i) = max_r;
  }

  double beta = std::max(arma::max(beta_min) + 1e-3, 1.0);

  // Update clusters using new beta
  List qt_x_result = qt_x_step_nystrom_cpp(n_rows, T, beta, B, col_sums, qy_t, qt);
  
  arma::mat qt_x_new = qt_x_result["qt_x"];
  bool rescue_occurred = qt_x_result["rescue_occurred"];

  return List::create(
    Named("qt_x") = qt_x_new,
    Named("beta") = beta,
    Named("rescue_occurred") = rescue_occurred
  );
}

// -----------------------------------------------------------------------------
// Nystrom version of qt_x_step_gib 
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat qt_x_step_gib_nystrom_cpp(int n_rows, int T, double beta, double alpha,
                                     const arma::mat& B,
                                     const arma::vec& col_sums,
                                     const arma::mat& qy_t,
                                     const arma::vec& qt) {

  // Initialize qt_x matrix 
  arma::mat qt_x(T, n_rows, arma::fill::zeros);

  // Precompute log(qt)
  arma::vec log_qt = arma::log(qt);
  log_qt.elem(arma::find(qt <= 1e-300)).fill(-1e100);

  // Precompute Cross-Entropy Matrix M
  arma::mat M = compute_cross_entropy_nystrom(B, col_sums, qy_t);

  for (int x = 0; x < n_rows; ++x) {
    // GIB score: log q(t) + beta * M(x,t)
    arma::vec scores = (log_qt + beta * M.row(x).t()) / alpha;
    
    // Numerical stability fix
    double max_score = scores.max();
    arma::vec l = arma::exp(scores - max_score);

    // Normalise
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
