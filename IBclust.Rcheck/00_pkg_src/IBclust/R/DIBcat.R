DIBcat <- function(X, ncl, randinit = NULL, lambda = -1,
                   maxiter = 100, nstart = 100,
                   verbose = FALSE) {

  # Validate inputs
  if (!is.data.frame(X)) {
    stop("Input 'X' must be a data frame.")
  }

  if (!is.numeric(ncl) || ncl <= 1 || ncl != round(ncl)) {
    stop("Input 'ncl' must be a positive integer greater than 1.")
  }

  if (!is.numeric(maxiter) || maxiter <= 0 || maxiter != round(maxiter)) {
    stop("'maxiter' must be a positive integer.")
  }

  if (!is.numeric(nstart) || nstart <= 0 || nstart != round(nstart)) {
    stop("'nstart' must be a positive integer.")
  }

  if (!is.null(randinit) && (!is.numeric(randinit) || length(randinit) != nrow(X))) {
    stop("'randinit' must be a numeric vector with length equal to the number of rows in 'X', or NULL.")
  }

  # Validate lambda
  if (!is.numeric(lambda) ||
      !(length(lambda) == 1 || length(lambda) == ncol(X)) ||
      any(lambda <= 0 & lambda != -1)) {
    stop("'lambda' must be either a single numeric value (-1 for automatic selection or a positive value) or a numeric vector with positive values matching the number of 'catcols'.")
  }

  # Additional check for maximum lambda value for nominal variables
  if (length(lambda) > 1 && length(lambda) == ncol(X)) {
    max_lambda <- sapply(1:ncol(X), function(col) {
      l <- length(unique(X[, col]))
      (l - 1) / l
    })
    if (any(lambda > max_lambda)) {
      stop("'lambda' values for nominal variables must not exceed their maximum allowable value of (l - 1)/l, where l is the number of categories in the variable.")
    }
  }

  # Preprocessing
  X <- preprocess_cat_data(X)

  if (length(lambda) == 1){
    if (lambda == -1){
      # Compute lambda for categorical data
      lambda <- compute_lambda_cat(X)
    }
  }

  bws_vec <- lambda
  
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
                               catcols = seq_len(ncol(X)), runs = nstart,
                               verbose = verbose)

  # Remove 's' element from the returned list
  best_clust$s <- NULL

  # Warning if clustering failed
  if (best_clust[[3]] == Inf) {
    warning("Initial cluster assignment remained unchanged; use other hyperparameter values for DIBcat to converge.")
  }

  return(best_clust)
}
