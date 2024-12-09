DIBcat <- function(X, ncl, randinit = NULL, lambda = -1,
                   maxiter = 100, nstart = 100, select_features = FALSE) {
  # Helper function to preprocess categorical data
  preprocess_cat_data <- function(X) {
    X <- data.frame(X)
    for (i in seq_len(ncol(X))) {
      X[, i] <- factor(X[, i], levels = unique(X[, i]), labels = seq_along(unique(X[, i])))
    }
    return(X)
  }

  # Helper function to compute lambda for categorical and ordinal data
  compute_lambda_cat <- function(X, lambda) {
    if (lambda == -1) {
      num_lvls_vec <- sapply(X, function(x) length(unique(x)))
      cat_rel_imp <- 2
      lambda <- (num_lvls_vec - 1) / (num_lvls_vec + cat_rel_imp - 1)
      ordvars <- sapply(c(1:ncol(X)), function(var) is.ordered(X[, var]))
      if (any(ordvars)) {
        inxs <- which(ordvars)
        lambda[inxs] <- (1 / cat_rel_imp)^(1 / (num_lvls_vec[inxs] - 1))
      }
    }
    return(lambda)
  }

  # Preprocessing
  X <- preprocess_cat_data(X)

  # Compute lambda for categorical data
  lambda <- compute_lambda_cat(X, lambda)

  # Feature selection (optional)
  if (select_features) {
    bw <- lambda  # Use lambda as bandwidth for categorical variables
    bws_vec <- eigengap(data = X, contcols = c(), catcols = seq_len(ncol(X)),
                        bw = bw, ncl = ncl)
  } else {
    bws_vec <- lambda
  }

  # Compute joint probability density for categorical variables
  pxy_list <- coord_to_pxy_R(as.data.frame(X), s = 0, cat_cols = seq_len(ncol(X)),
                             cont_cols = c(), lambda = bws_vec)
  py_x <- pxy_list$py_x
  px <- pxy_list$px
  hy <- pxy_list$hy

  # Run DIB iteration for clustering
  best_clust <- DIBmix_iterate(X, ncl = ncl, randinit = randinit, tol = 0,
                               py_x = py_x, hy = hy, px = px, maxiter = maxiter,
                               bws_vec = bws_vec, contcols = c(),
                               catcols = seq_len(ncol(X)), runs = nstart)

  # Remove 's' element from the returned list
  best_clust$s <- NULL

  # Warning if clustering failed
  if (best_clust[[3]] == Inf) {
    warning("Initial cluster assignment remained unchanged; use other hyperparameter values for DIBcat to converge.")
  }

  return(best_clust)
}
