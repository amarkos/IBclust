DIBmix <- function(X, ncl, catcols, contcols, randinit = NULL,
                   lambda = -1, s = -1, scale = TRUE,
                   maxiter = 100, nstart = 100,
                   select_features = FALSE) {

  # Validate inputs
  if (!is.data.frame(X)) {
    stop("Input 'X' must be a data frame.")
  }

  if (!is.numeric(ncl) || ncl <= 1 || ncl != round(ncl)) {
    stop("Input 'ncl' must be a positive integer greater than 1.")
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

  if (!is.logical(select_features)) {
    stop("'select_features' must be a logical value (TRUE or FALSE).")
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
  # Helper function to preprocess data
  preprocess_data <- function(X, catcols, contcols) {
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
    return(X)
  }

  # Helper function to compute s using the bandwidth heuristic
  compute_bandwidth <- function(X, contcols, s) {
    if (s == -1) {
      s_seq <- seq(0.1, 10, by = 1e-1)
      for (s_val in s_seq) {
        pxy_list_cont <- coord_to_pxy_R(as.data.frame(X[, contcols]), s = s_val,
                                        cat_cols = c(), cont_cols = seq_along(contcols),
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

  # Helper function to compute lambda
  compute_lambda <- function(X, catcols, contcols, avg_max_min_py_x, cat_ratio, cat_rel_imp) {
    num_lvls_vec <- sapply(X[, catcols, drop = FALSE], function(x) length(unique(x)))
    cat_rel_imp <- min(2, avg_max_min_py_x^(1 / length(contcols)))
    if (cat_rel_imp == 2 & cat_ratio > 0.5) {
      cat_rel_imp <- 2 - cat_ratio
    }
    lambda <- (num_lvls_vec - 1) / (num_lvls_vec + cat_rel_imp - 1)
    ordvars <- sapply(catcols, function(var) is.ordered(X[, var]))
    if (any(ordvars)) {
      inxs <- which(ordvars)
      lambda[inxs] <- (1 / cat_rel_imp)^(1 / (num_lvls_vec[inxs] - 1))
    }
    return(lambda)
  }

  # Preprocessing
  cat_ratio <- length(catcols) / ncol(X)
  cat_first <- FALSE
  X <- preprocess_data(X, catcols, contcols)

  # Bandwidth computation
  if (length(s) == 1)
    if (s == -1)
      s <- compute_bandwidth(X, contcols, s)


  # Compute py_x for continuous variables
  pxy_list_cont <- coord_to_pxy_R(as.data.frame(X[, contcols]), s = s, cat_cols = c(),
                                  cont_cols = seq_along(contcols), lambda = 0)
  pyx_cont <- pxy_list_cont$py_x
  avg_max_min_py_x <- mean(apply(pyx_cont, 2, function(x) max(x) / min(x)))

  if (length(lambda) == 1)
    if (lambda == -1)
      # Compute lambda
      lambda <- compute_lambda(X, catcols, contcols, avg_max_min_py_x, cat_ratio, cat_rel_imp)

  # Final bandwidth values
  bw <- numeric(ncol(X))
  bw[contcols] <- s
  bw[catcols] <- lambda

  # Feature selection
  if (select_features) {
    bws_vec <- eigengap(data = X, contcols = contcols, catcols = catcols, bw = bw, ncl = ncl)
  } else {
    bws_vec <- bw
  }
  print(bws_vec[catcols])
  # Construct joint density with final bandwidths
  pxy_list <- coord_to_pxy_R(X, s = bws_vec[contcols], cat_cols = catcols,
                             cont_cols = contcols, lambda = bws_vec[catcols])
  py_x <- pxy_list$py_x
  px <- pxy_list$px
  hy <- pxy_list$hy

  # Run DIBmix iteration
  best_clust <- DIBmix_iterate(X, ncl = ncl, randinit = randinit, tol = 0,
                               py_x = py_x, hy = hy, px = px, maxiter = maxiter,
                               bws_vec = bws_vec, contcols = contcols,
                               catcols = catcols, runs = nstart)

  # Warning if clustering failed
  if (best_clust$MutualInfo == Inf) {
    warning("Initial cluster assignment remained unchanged; use other hyperparameter values for DIBmix to converge.")
  }

  return(best_clust)
}
