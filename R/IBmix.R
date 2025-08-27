IBmix <- function(X, ncl, beta, randinit = NULL,
                  s = -1, lambda = -1, scale = TRUE,
                  maxiter = 100, nstart = 100, contkernel = "gaussian",
                  nomkernel = "aitchisonaitken", ordkernel = "liracine",
                  cat_first = FALSE, verbose = FALSE) {
  
  # Validate inputs
  if (!is.data.frame(X)) {
    stop("Input 'X' must be a data frame.")
  }
  if (!is.numeric(ncl) || ncl <= 1 || ncl != round(ncl)) {
    stop("Input 'ncl' must be a positive integer greater than 1.")
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
  if (!is.logical(cat_first)) {
    stop("'cat_first' must be a logical value (TRUE or FALSE).")
  }
  # Check kernel types
  if (!contkernel %in% c("gaussian", "epanechnikov")){
    stop("'contkernel' can only be one of 'gaussian' or 'epanechnikov'")
  }
  if (!nomkernel %in% c("aitchisonaitken", "liracine")){
    stop("'nomkernel' can only be one of 'aitchisonaitken' or 'liracine'")
  }
  if (!ordkernel %in% c("liracine", "wangvanryzin")){
    stop("'ordkernel' can only be one of 'liracine' or 'wangvanryzin'")
  }
  
  X <- data.frame(X)
  # Check catcols/contcols
  catcols <- as.integer(which(sapply(X, is.factor)))
  if (length(catcols) == 0){
    catcols <- c()
  }
  contcols <- as.integer(which(sapply(X, is.numeric)))
  if (length(contcols) == 0){
    contcols <- c()
  }
  
  if (cat_first & any(c(s, lambda) != -1)){
    stop("'cat_first' can only be TRUE when all bandwidths are determined by the algorithm (s = -1, lambda = -1).")
  }
  if (cat_first & length(catcols) == 0){
    stop("'cat_first' can only be TRUE when there are categorical variables in the data set.")
  }
  
  # Validate lambda if any categorical columns exist
  if (length(catcols) > 0){
    if (!is.numeric(lambda) ||
        !(length(lambda) == 1 || length(lambda) == length(catcols)) ||
        any(lambda <= 0 & lambda != -1)) {
      stop("'lambda' must be either a single numeric value (-1 for automatic selection or a positive value) or a numeric vector with positive values matching the number of 'catcols'.")
    }
    # Additional check for maximum lambda value for nominal variables
    if (length(lambda) > 1 && length(lambda) == length(catcols)) {
      if (nomkernel == "liracine"){
        max_lambda <- sapply(catcols, function(col) {
          l <- length(unique(X[, col]))
          (l - 1) / l
        })
      } else {
        max_lambda <- 1
      }
      if (any(lambda > max_lambda)) {
        stop("'lambda' values for nominal variables must not exceed their maximum allowable value.")
      }
    }
    X[, catcols] <- preprocess_cat_data(X[, catcols])
  }
  
  # Validate s
  if (length(contcols) > 0){
    if (!is.numeric(s) ||
        !(length(s) == 1 || length(s) == length(contcols)) ||
        any(s <= 0 & s != -1)) {
      stop("'s' must be either a single numeric value (-1 for automatic selection or a positive value) or a numeric vector with positive values matching the number of 'contcols'.")
    }
    if (scale){
      X[, contcols] <- as.data.frame(preprocess_cont_data(X[, contcols]))
    }
  }
  
  if (length(contcols) == 0){
    if (length(lambda) == 1){
      if (lambda == -1){
        # Compute lambda for categorical data
        lambda <- compute_lambda_cat(X, nomkernel, ordkernel)
      }
    }
    bws_vec <- lambda
  } else if (length(catcols) == 0){
    if (length(s) == 1){
      if (s == -1){
        s <- compute_bandwidth_cont(X, contkernel = contkernel)
      }
    }
    bws_vec <- rep(s, length(contcols))
  } else {
    bws_vec <- compute_s_lambda(X, contcols, catcols, s, lambda,
                                contkernel, nomkernel, ordkernel,
                                cat_first)
  }
  
  # Construct joint density with final bandwidths
  # Construct joint density with final bandwidths
  pxy_list <- coord_to_pxy_R(as.data.frame(X),
                             s = if (length(contcols) > 0){
                               bws_vec[contcols]
                             } else {
                               -1
                             },
                             lambda = if (length(catcols) > 0){
                               bws_vec[catcols]
                             } else {
                               -1
                             },
                             cat_cols = catcols,
                             cont_cols = contcols,
                             contkernel = contkernel,
                             nomkernel = nomkernel,
                             ordkernel = ordkernel)
  
  py_x <- pxy_list$py_x
  px <- pxy_list$px
  pxy <- pxy_list$pxy
  hy <- pxy_list$hy
  
  ######################################################
  best_clust <- IBmix_iterate(X, ncl = ncl, beta = beta,
                              randinit = randinit,
                              tol = 0, py_x, hy, px, maxiter,
                              bws_vec, contcols, catcols,
                              runs = nstart, verbose = verbose)
  ######################################################
  
  return(best_clust)
}
