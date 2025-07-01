# Helper function to preprocess continuous data
preprocess_cont_data <- function(X) {
  X <- data.frame(X)
  X <- scale(X)  # Standardize continuous variables
  return(X)
}

# Helper function to preprocess categorical data
preprocess_cat_data <- function(X) {
  X <- data.frame(X)
  for (i in seq_len(ncol(X))) {
    X[, i] <- factor(X[, i], levels = unique(X[, i]), labels = seq_along(unique(X[, i])))
  }
  return(X)
}

# Helper function to compute bandwidth (s) for continuous data
compute_bandwidth_cont <- function(X){
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
  return(s_val)
}

# Helper function to compute lambda for categorical and ordinal data
compute_lambda_cat <- function(X){
  num_lvls_vec <- sapply(X, function(x) length(unique(x)))
  cat_rel_imp <- 2
  lambda <- (num_lvls_vec - 1) / (num_lvls_vec + cat_rel_imp - 1)
  ordvars <- sapply(c(1:ncol(X)), function(var) is.ordered(X[, var]))
  if (any(ordvars)) {
    inxs <- which(ordvars)
    lambda[inxs] <- (1 / cat_rel_imp)^(1 / (num_lvls_vec[inxs] - 1))
  }
  return(lambda)
}

# Helper function to compute s and lambda to equalise variable contributions
compute_s_lambda <- function(X, contcols, catcols, s, lambda){
  cat_rel_imp <- 2
  # Get ratio of categorical variables
  cat_ratio <- length(catcols)/(ncol(X))
  # Indicator for whether \zeta_j is determined by \ksi_j
  cat_first <- FALSE
  
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
  bw[catcols] <- lambda
  bws_vec <- bw
  return(bws_vec)
}
