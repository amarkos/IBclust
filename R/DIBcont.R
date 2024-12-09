DIBcont <- function(X, ncl, randinit = NULL, s = -1, scale = TRUE,
                    maxiter = 100, nstart = 100, select_features = FALSE) {
  # Helper function to preprocess data
  preprocess_cont_data <- function(X) {
    X <- data.frame(X)
    X <- scale(X)  # Standardize continuous variables
    return(X)
  }

  # Helper function to compute bandwidth (s) for continuous data
  compute_bandwidth_cont <- function(X, s) {
    if (s == -1) {
      s_seq <- seq(0.1, 10, by = 1e-1)
      for (s_val in s_seq) {
        pxy_list_cont <- coord_to_pxy_R(as.data.frame(X), s = s_val,
                                        cat_cols = c(), cont_cols = seq_len(ncol(X)),
                                        lambda = 0)
        pyx_cont <- pxy_list_cont$py_x
        avg_py_x <- mean(apply(pyx_cont, 2, function(x) max(x) / max(x[-which.max(x)])))
        if (avg_py_x < 1.1) {
          return(s_val - 1e-1)
        }
      }
    }
    return(s)
  }

  # Preprocessing
  if (scale == FALSE)
    X <- preprocess_cont_data(X)

  # Bandwidth computation
  s <- compute_bandwidth_cont(X, s)

  # Compute joint probability density for continuous variables
  pxy_list <- coord_to_pxy_R(as.data.frame(X), s = s, cat_cols = c(),
                             cont_cols = seq_len(ncol(X)), lambda = 0)
  py_x <- pxy_list$py_x
  px <- pxy_list$px
  hy <- pxy_list$hy

  # Feature selection using eigengap heuristic (optional)
  if (select_features) {
    bw <- rep(s, ncol(X))
    bws_vec <- eigengap(data = X, contcols = seq_len(ncol(X)), catcols = c(),
                        bw = bw, ncl = ncl)
  } else {
    bws_vec <- rep(s, ncol(X))
  }

  # Run DIB iteration for clustering
  best_clust <- DIBmix_iterate(X, ncl = ncl, randinit = randinit, tol = 0,
                               py_x = py_x, hy = hy, px = px, maxiter = maxiter,
                               bws_vec = bws_vec, contcols = seq_len(ncol(X)),
                               catcols = c(), runs = nstart)

  # Warning if clustering failed
  if (best_clust[[3]] == Inf) {
    warning("Initial cluster assignment remained unchanged; use other hyperparameter values for DIBcont to converge.")
  }

  return(best_clust)
}
