IBcat <- function(X, ncl, beta, randinit = NULL, lambda = -1,
                  maxiter = 100, nstart = 100, select_features = FALSE,
                  verbose = FALSE) {
  
  # Validate inputs
  if (!is.data.frame(X)) {
    stop("Input 'X' must be a data frame.")
  }
  
  if (!is.numeric(ncl) || ncl <= 1 || ncl != round(ncl)) {
    stop("Input 'ncl' must be a positive integer greater than 1.")
  }
  
  if (!is.numeric(beta) || beta <= 0) {
    stop("Input 'beta' must be a positive number.")
  }
  
  if (!is.numeric(maxiter) || maxiter <= 0 || maxiter != round(maxiter)) {
    stop("'maxiter' must be a positive integer.")
  }
  
  if (!is.numeric(nstart) || nstart <= 0 || nstart != round(nstart)) {
    stop("'nstart' must be a positive integer.")
  }
  
  if (!is.logical(select_features)) {
    stop("'select_features' must be a logical value (TRUE or FALSE).")
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
  
  if (length(lambda) == 1)
    if (lambda == -1)
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
  best_clust <- IBmix_iterate(X, ncl = ncl, beta = beta, randinit = randinit, tol = 0,
                              py_x = py_x, hy = hy, px = px, maxiter = maxiter,
                              bws_vec = bws_vec, contcols = c(),
                              catcols = seq_len(ncol(X)), runs = nstart,
                              verbose = verbose)
  
  # Remove 's' element from the returned list
  best_clust$s <- NULL
  
  return(best_clust)
}
