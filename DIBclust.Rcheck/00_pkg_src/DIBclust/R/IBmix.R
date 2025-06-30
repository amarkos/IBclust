IBmix <- function(X, ncl, beta, catcols, contcols, randinit = NULL,
                  lambda = -1, s = -1, scale = TRUE,
                  maxiter = 100, nstart = 100,
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
  
  if (!all(catcols %in% seq_along(X))) {
    stop("Some 'catcols' indices are out of bounds or invalid.")
  }
  
  if (!all(contcols %in% seq_along(X))) {
    stop("Some 'contcols' indices are out of bounds or invalid.")
  }
  
  if (any(duplicated(c(catcols, contcols)))) {
    stop("'catcols' and 'contcols' must not overlap.")
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
  
  # Validate lambda
  if (!is.numeric(lambda) ||
      !(length(lambda) == 1 || length(lambda) == length(catcols)) ||
      any(lambda <= 0 & lambda != -1)) {
    stop("'lambda' must be either a single numeric value (-1 for automatic selection or a positive value) or a numeric vector with positive values matching the number of 'catcols'.")
  }
  
  # Additional check for maximum lambda value for nominal variables
  if (length(lambda) > 1 && length(lambda) == length(catcols)) {
    max_lambda <- sapply(catcols, function(col) {
      l <- length(unique(X[, col]))
      (l - 1) / l
    })
    if (any(lambda > max_lambda)) {
      stop("'lambda' values for nominal variables must not exceed their maximum allowable value of (l - 1)/l, where l is the number of categories in the variable.")
    }
  }
  
  
  # Validate s
  if (!is.numeric(s) ||
      !(length(s) == 1 || length(s) == length(contcols)) ||
      any(s <= 0 & s != -1)) {
    stop("'s' must be either a single numeric value (-1 for automatic selection or a positive value) or a numeric vector with positive values matching the number of 'contcols'.")
  }
  
  cat_rel_imp <- 2
  # Get ratio of categorical variables
  cat_ratio <- length(catcols)/(ncol(X))
  # Indicator for whether \zeta_j is determined by \ksi_j
  cat_first <- FALSE
  
  X <- data.frame(X)
  for (i in catcols) {
    X[, i] <- factor(X[, i], levels = unique(X[, i]), labels = seq_along(unique(X[, i])))
  }
  if (length(contcols) == 1) {
    if (scale == TRUE)
      X[, contcols] <- as.numeric(scale(X[, contcols]))
  } else {
    if (scale == TRUE)
      X[, contcols] <- scale(X[, contcols])
  }
  
  # Bandwidth values
  if (length(s) == 1){
    if (s == -1){
      s_seq <- seq(0.1, 10, by=1e-1)
      for (s_val in s_seq){
        pxy_list_cont <- coord_to_pxy_R(as.data.frame(X[, contcols]), s = s_val, cat_cols = c(),
                                        cont_cols = c(1:length(contcols)),
                                        lambda = 0)
        pyx_cont <- pxy_list_cont$py_x
        avg_py_x <- mean(apply(pyx_cont, 2, FUN = function(x) max(x)/max(x[-which.max(x)])))
        if (avg_py_x < (1.1)){
          s <- s_val - 1e-1
          avg_max_min_py_x <- mean(apply(pyx_cont, 2, FUN = function(x) max(x)/min(x)))
          break
        }
      }
    }
  } else {
    pxy_list_cont <- coord_to_pxy_R(as.data.frame(X[, contcols]), s = s, cat_cols = c(),
                                    cont_cols = c(1:length(contcols)),
                                    lambda = 0)
    pyx_cont <- pxy_list_cont$py_x
  }
  if (length(lambda) == 1 && lambda == -1){
    num_lvls_vec <- sapply(X[, catcols, drop = FALSE], function(x) length(unique(x)))
    cat_rel_imp <- min(2, avg_max_min_py_x^(1/length(contcols)))
    # Check if cat_rel_imp is 2 & cat_ratio > 0.5
    if (cat_rel_imp == 2 & cat_ratio > 0.5){
      cat_rel_imp <- 2 - cat_ratio
      cat_first <- TRUE
    }
    ##lambda <- (num_lvls_vec-1)/(num_lvls_vec)
    lambda <- (num_lvls_vec - 1)/(num_lvls_vec + cat_rel_imp - 1)
    # Check for ordinal variables
    ordvars <- c()
    for (var in catcols){
      if (is.ordered(X[, var])){
        ordvars <- c(ordvars, var)
      }
    }
    if (length(ordvars) > 0){
      inxs <- match(ordvars, catcols)
      lambda[inxs] <- (1/cat_rel_imp)^(1/(num_lvls_vec[inxs]-1))
    }
  }
  
  if (cat_first){
    s_seq <- seq(0.1, 10, by=1e-1)
    for (s_val in s_seq){
      pxy_list_cont <- coord_to_pxy_R(as.data.frame(X[, contcols]), s = s_val, cat_cols = c(),
                                      cont_cols = c(1:length(contcols)),
                                      lambda = 0)
      pyx_cont <- pxy_list_cont$py_x
      avg_max_min_py_x <- mean(apply(pyx_cont, 2, FUN = function(x) max(x)/min(x)))
      if (avg_max_min_py_x < (cat_rel_imp)^(length(catcols))){
        s <- s_val
        break
      }
    }
  }
  
  bw <- rep(NA, ncol(X))
  bw[contcols] <- s
  ##bw[catcols] <- lambda - eps_star
  bw[catcols] <- lambda
  bws_vec <- bw
  
  # Construct joint density with final bandwidths
  pxy_list <- coord_to_pxy_R(X, s = bws_vec[contcols],
                             cat_cols = catcols, cont_cols = contcols,
                             lambda = bws_vec[catcols])
  
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
