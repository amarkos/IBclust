DIBmix <- function(X, ncl, catcols, contcols, randinit = NULL,
                   lambda = -1, s = -1, scale = TRUE,
                   maxiter = 100, nstart = 100,
                   select_features = FALSE) {

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
  s <- compute_bandwidth(X, contcols, s)

  # Compute py_x for continuous variables
  pxy_list_cont <- coord_to_pxy_R(as.data.frame(X[, contcols]), s = s, cat_cols = c(),
                                  cont_cols = seq_along(contcols), lambda = 0)
  pyx_cont <- pxy_list_cont$py_x
  avg_max_min_py_x <- mean(apply(pyx_cont, 2, function(x) max(x) / min(x)))

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
