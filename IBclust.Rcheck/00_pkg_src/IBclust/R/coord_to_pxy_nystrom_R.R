#' Compute pxy using Nyström Approximation
#'
#' This function computes an approximate pxy matrix using Nyström method with landmark points.
#'
#' @param X Data frame containing the data set.
#' @param s Bandwidth parameter for continuous variables.
#' @param cat_cols Indices of categorical columns in X.
#' @param cont_cols Indices of continuous columns in X.
#' @param lambda Bandwidth parameter for categorical variables.
#' @param n_landmarks Number of landmark points (default: ceiling(sqrt(nrow(X)))).
#' @param contkernel Continuous kernel (Gaussian or Epanechnikov)
#' @param nomkernel Unordered categorical (nominal) kernel (Aitchison & Aitken or Li & Racine)
#' @param ordkernel Ordered categorical (ordinal) kernel (Li & Racine or Wang & van Ryzin)
#'
#' @return A list containing py_x_approx (Nyström approximation), px, pxy, hy, and landmark_indices.
#'
#' @keywords internal
#' @noRd
coord_to_pxy_nystrom_R <- function(X, s, cat_cols, cont_cols, lambda,
                                   n_landmarks = NULL,
                                   contkernel = "gaussian",
                                   nomkernel = "aitchisonaitken",
                                   ordkernel = "liracine"){
  
  n <- nrow(X)
  # Default number of landmarks
  if (is.null(n_landmarks)) {
    n_landmarks <- min(ceiling(sqrt(n)), n %/% 2)
  }
  # Ensure valid number of landmarks
  if (n_landmarks >= n) {
    warning("n_landmarks >= n, using full matrix computation instead")
    return(coord_to_pxy_R(X, s, cat_cols, cont_cols, lambda,
                          contkernel, nomkernel, ordkernel))
  }
  
  # Select random landmark points
  landmark_idx <- sample(n, n_landmarks, replace = FALSE)
  X_landmarks <- X[landmark_idx, , drop = FALSE]
  
  # Prepare bandwidth vector
  bws <- rep(NA, ncol(X))
  bws[cat_cols] <- lambda
  bws[cont_cols] <- s
  
  # C: k×k kernel matrix between landmarks
  C <- t(np::npksum(bws = bws,
                    txdat = X_landmarks,
                    exdat = X_landmarks,
                    ckertype = contkernel,
                    ukertype = nomkernel,
                    okertype = ordkernel,
                    return.kernel.weights = TRUE)$kw)
  
  # W: n×k kernel matrix (all points × landmarks)
  W <- t(np::npksum(bws = bws,
                    txdat = X_landmarks,
                    exdat = X,
                    ckertype = contkernel,
                    ukertype = nomkernel,
                    okertype = ordkernel,
                    return.kernel.weights = TRUE)$kw)
  
  C <- sweep(C, 2, colSums(C), '/')
  W <- sweep(W, 2, colSums(W), '/')
  
  # Compute inverse sqrt of C using svd
  C_svd <- svd(C)
  tol <- max(dim(C)) * .Machine$double.eps * max(C_svd$d)
  pos <- C_svd$d > tol
  
  if (sum(pos) < length(C_svd$d)) {
    warning(sprintf("C matrix is rank-deficient: rank %d vs. dimension %d",
                    sum(pos), length(C_svd$d)))
  }
  
  C_sqrt_inv <- C_svd$v[, pos, drop = FALSE] %*% 
    diag(1/sqrt(C_svd$d[pos]), nrow = sum(pos)) %*%
    t(C_svd$u[, pos, drop = FALSE])
  
  # B factor: n × r matrix
  B <- W %*% C_sqrt_inv
  
  col_sums_B <- colSums(B)
  col_sums_K <- as.vector(B %*% col_sums_B)
  col_sums_K <- pmax(col_sums_K, 1e-8)
  
  # Check for problematic columns
  if (any(col_sums_K <= 0)) {
    warning(sprintf("Found %d columns with non-positive sums", 
                    sum(col_sums_K <= 0)))
    col_sums_K[col_sums_K <= 0] <- 1e-10
  }
  # Create factored representation object
  py_x_nystrom <- list(
    B = B,
    col_sums = col_sums_K,
    n = n,
    r = ncol(B),
    landmark_indices = landmark_idx
  )
  class(py_x_nystrom) <- "nystrom_matrix"
  
  D_vec <- 1 / col_sums_K
  row_sums_py_x <- as.vector(B %*% (t(B) %*% D_vec))
  pxy_col_sums <- row_sums_py_x / n
  pxy_col_sums <- pmax(pxy_col_sums, 0)
  
  # Renormalise to sum to 1
  pxy_col_sums <- pxy_col_sums / sum(pxy_col_sums)
  
  # Compute entropy
  hy <- entropySingle(pxy_col_sums)
  px <- matrix(1/n, nrow = n, ncol = 1)
  
  return(list('py_x' = py_x_nystrom,
              'px' = px,
              'pxy' = NULL,
              'hy' = hy,
              'landmark_indices' = landmark_idx,
              'n_landmarks' = n_landmarks,
              'is_nystrom' = TRUE))
}
