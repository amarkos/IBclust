GIBcont <- function(X, ncl, beta, alpha, randinit = NULL, s = -1, scale = TRUE,
                    maxiter = 100, nstart = 100,
                    verbose = FALSE) {
  
  # Validate inputs
  if (!is.numeric(ncl) || ncl <= 1 || ncl != round(ncl)) {
    stop("Input 'ncl' must be a positive integer greater than 1.")
  }
  
  if (!is.numeric(beta) || beta <= 0) {
    stop("Input 'beta' must be a positive number.")
  }
  
  if (!is.numeric(alpha) || alpha < 0) {
    stop("Input 'alpha' must be a non-negative number.")
  }
  
  if (!is.logical(scale)) {
    stop("'scale' must be a logical value (TRUE or FALSE).")
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
  
  # Validate s
  if (!is.numeric(s) ||
      !(length(s) == 1 || length(s) == ncol(X)) ||
      any(s <= 0 & s != -1)) {
    stop("'s' must be either a single numeric value (-1 for automatic selection or a positive value) or a numeric vector with positive values matching the number of 'contcols'.")
  }
  
  # Check special case of alpha = 0 (DIBmix) or alpha = 1 (IBmix)
  if (alpha == 1){
    message('alpha = 1; running IBcont.')
    best_clust <- IBcont(X, ncl, beta, randinit,
                        s, scale, maxiter, nstart,
                        verbose)
  } else if (alpha == 0){
    message('alpha = 0; running DIBcont - value of beta is ignored.')
    best_clust <- DIBcont(X, ncl, randinit,
                         s, scale, maxiter, nstart,
                         verbose)
  } else {
    # Preprocessing
    if (scale){
      X <- preprocess_cont_data(X)
    }
    
    # Bandwidth computation
    if (length(s) == 1){
      if (s == -1){
        s <- compute_bandwidth_cont(X)
      }
    }
    
    # Compute joint probability density for continuous variables
    pxy_list <- coord_to_pxy_R(as.data.frame(X), s = s, cat_cols = c(),
                               cont_cols = seq_len(ncol(X)), lambda = 0)
    py_x <- pxy_list$py_x
    px <- pxy_list$px
    hy <- pxy_list$hy
    
    bws_vec <- rep(s, ncol(X))
    
    # Run GIB iteration for clustering
    best_clust <- GIBmix_iterate(X, ncl = ncl, beta = beta, alpha = alpha, randinit = randinit, tol = 0,
                                 py_x = py_x, hy = hy, px = px, maxiter = maxiter,
                                 bws_vec = bws_vec, contcols = seq_len(ncol(X)),
                                 catcols = c(), runs = nstart, verbose = verbose)
  }
  return(best_clust)
}
