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
compute_bandwidth_cont <- function(X, contkernel){
  s_seq <- seq(0.1, 10, by = 1e-1)
  for (s_val in s_seq) {
    pxy_list_cont <- coord_to_pxy_R(as.data.frame(X), s = s_val,
                                    cat_cols = c(), cont_cols = seq_len(ncol(X)),
                                    lambda = 0,
                                    contkernel = contkernel)
    pyx_cont <- pxy_list_cont$py_x
    nearest_neighbours_ratios <- apply(pyx_cont, 2, FUN = function(x) max(x)/max(x[-which.max(x)]))
    if (any(nearest_neighbours_ratios == Inf)) next
    avg_py_x <- mean(nearest_neighbours_ratios)
    ratio_thresh <- ifelse(contkernel == "gaussian", 1.1, 1.04)
    if (avg_py_x < ratio_thresh) {
      return(s_val - 1e-1)
    }
  }
  return(s_val)
}

# Helper function to compute lambda for categorical and ordinal data
compute_lambda_cat <- function(X, nomkernel, ordkernel){
  num_lvls_vec <- sapply(X, function(x) length(unique(x)))
  cat_rel_imp <- 2
  if (nomkernel == "aitchisonaitken"){
    lambda <- (num_lvls_vec - 1) / (num_lvls_vec + cat_rel_imp - 1)
  } else if (nomkernel == "liracine"){
    lambda <- rep(1/cat_rel_imp, length(num_lvls_vec))
  }
  
  ordvars <- sapply(c(1:ncol(X)), function(var) is.ordered(X[, var]))
  if (any(ordvars)) {
    inxs <- which(ordvars)
    if (ordkernel == "liracine"){
      lambda[inxs] <- (1 / cat_rel_imp)^(1 / (num_lvls_vec[inxs] - 1))
    } else if (ordkernel == "wangvanryzin"){
      cat_rel_imp_wvr <- 2*cat_rel_imp
      lambda[inxs] <- (1 / (cat_rel_imp_wvr/2))^(1 / (num_lvls_vec[inxs] - 1))
    }
    
  }
  return(lambda)
}

# Helper function to compute s and lambda to equalise variable contributions
compute_s_lambda <- function(X, contcols, catcols, s, lambda,
                             contkernel, nomkernel, ordkernel,
                             cat_first){
  cat_rel_imp <- 2
  # Get ratio of categorical variables
  cat_ratio <- length(catcols)/(ncol(X))
  
  # Check for ordinal variables
  ordvars <- c()
  for (var in catcols){
    if (is.ordered(X[, var])){
      ordvars <- c(ordvars, var)
    }
  }
  num_lvls_vec <- sapply(X[, catcols, drop = FALSE], function(x) length(unique(x)))
  
  # Bandwidth values
  if (!cat_first){
    if (length(s) == 1){
      if (s == -1){
        s_seq <- seq(0.1, 10, by=1e-1)
        for (s_val in s_seq){
          pxy_list_cont <- coord_to_pxy_R(as.data.frame(X[, contcols]), s = s_val, cat_cols = c(),
                                          cont_cols = c(1:length(contcols)),
                                          lambda = 0,
                                          contkernel = contkernel,
                                          nomkernel = nomkernel,
                                          ordkernel = ordkernel)
          pyx_cont <- pxy_list_cont$py_x
          # Filter out infinites caused by division by zero (e.g. Epanechnikov kernel)
          nearest_neighbours_ratios <- apply(pyx_cont, 2, FUN = function(x) max(x)/max(x[-which.max(x)]))
          if (any(nearest_neighbours_ratios == Inf)) next
          avg_py_x <- mean(nearest_neighbours_ratios)
          ratio_thresh <- ifelse(contkernel == "gaussian", 1.1, 1.04)
          if (avg_py_x < ratio_thresh){
            s <- s_val - 1e-1
            avg_max_min_py_x <- mean(apply(pyx_cont, 2, FUN = function(x) max(x)/min(x[which(x > 0)])))
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
      cat_rel_imp <- min(2, avg_max_min_py_x^(1/length(contcols)))
      # Check if cat_rel_imp is 2 & cat_ratio > 0.5
      if (cat_rel_imp == 2 & cat_ratio > 0.5){
        cat_rel_imp <- 2 - cat_ratio
        cat_first <- TRUE
      }
      if (nomkernel == "aitchisonaitken"){
        lambda <- (num_lvls_vec - 1) / (num_lvls_vec + cat_rel_imp - 1)
      } else if (nomkernel == "liracine"){
        lambda <- rep(1/cat_rel_imp, length(num_lvls_vec))
      }
      if (length(ordvars) > 0){
        inxs <- match(ordvars, catcols)
        if (ordkernel == "liracine"){
          lambda[inxs] <- (1 / cat_rel_imp)^(1 / (num_lvls_vec[inxs] - 1))
        } else if (ordkernel == "wangvanryzin"){
          cat_rel_imp_wvr <- 2*(2 - cat_ratio)
          lambda[inxs] <- (1 / (cat_rel_imp_wvr/2))^(1 / (num_lvls_vec[inxs] - 1))
        }
      }
    }
  }
  if (cat_first){
    cat_rel_imp <- 2 - cat_ratio
    cat_first <- TRUE
    if (nomkernel == "aitchisonaitken"){
      lambda <- (num_lvls_vec - 1) / (num_lvls_vec + cat_rel_imp - 1)
    } else if (nomkernel == "liracine"){
      lambda <- rep(1/cat_rel_imp, length(num_lvls_vec))
    }
    if (length(ordvars) > 0){
      inxs <- match(ordvars, catcols)
      if (ordkernel == "liracine"){
        lambda[inxs] <- (1 / cat_rel_imp)^(1 / (num_lvls_vec[inxs] - 1))
      } else if (ordkernel == "wangvanryzin"){
        cat_rel_imp_wvr <- 2*(2 - cat_ratio)
        lambda[inxs] <- (1 / (cat_rel_imp_wvr/2))^(1 / (num_lvls_vec[inxs] - 1))
      }
    }
    s_seq <- seq(0.1, 10, by=1e-1)
    for (s_val in s_seq){
      pxy_list_cont <- coord_to_pxy_R(as.data.frame(X[, contcols]), s = s_val, cat_cols = c(),
                                      cont_cols = c(1:length(contcols)),
                                      lambda = 0, contkernel, nomkernel, ordkernel)
      pyx_cont <- pxy_list_cont$py_x
      avg_max_min_py_x <- mean(apply(pyx_cont, 2, FUN = function(x) max(x)/min(x[which(x > 0)])))
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
