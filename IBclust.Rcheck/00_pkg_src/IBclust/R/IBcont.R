IBcont <- function(X, ncl, beta, randinit = NULL, s = -1, scale = TRUE,
                   maxiter = 100, nstart = 100,
                   verbose = FALSE) {
  
  # Validate inputs
  if (!is.numeric(ncl) || ncl <= 1 || ncl != round(ncl)) {
    stop("Input 'ncl' must be a positive integer greater than 1.")
  }
  
  if (!is.numeric(beta) || beta <= 0) {
    stop("Input 'beta' must be a positive number.")
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
  
  # Run IB iteration for clustering
  best_clust <- IBmix_iterate(X, ncl = ncl, beta = beta, randinit = randinit, tol = 0,
                              py_x = py_x, hy = hy, px = px, maxiter = maxiter,
                              bws_vec = bws_vec, contcols = seq_len(ncol(X)),
                              catcols = c(), runs = nstart, verbose = verbose)
  
  return(best_clust)
}
