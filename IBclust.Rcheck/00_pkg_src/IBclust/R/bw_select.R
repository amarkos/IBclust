find_s_nn <- function(X, contcols, contkernel, nomkernel, ordkernel,
                      n_bw_reps = 100, s_min = 0.1, s_max = 100,
                      nystrom) {
  ratio_thresh <- ifelse(contkernel == "gaussian", 1.1, 1.04)
  if (nystrom) {
    n_bw_sample <- min(nrow(X), 1000)
    compute_avg <- function(s_val) {
      avg_py_x_reps <- replicate(n_bw_reps, {
        bw_sample_idx <- sample(nrow(X), n_bw_sample)
        X_bw <- X[bw_sample_idx, ]
        
        pxy_list_cont <- coord_to_pxy_R(
          as.data.frame(X_bw[, contcols]),
          s = s_val,
          cat_cols = c(),
          cont_cols = seq_along(contcols),
          lambda = 0,
          contkernel = contkernel,
          nomkernel = nomkernel,
          ordkernel = ordkernel
        )
        pyx_cont <- pxy_list_cont$py_x
        nearest_neighbours_ratios <- apply(pyx_cont, 2, function(x) max(x) / max(x[-which.max(x)]))
        mean(nearest_neighbours_ratios)
      })
      median(avg_py_x_reps[is.finite(avg_py_x_reps)])
    }
  } else {
    compute_avg <- function(s_val) {
      pxy_list_cont <- coord_to_pxy_R(
        as.data.frame(X[, contcols]),
        s = s_val,
        cat_cols = c(),
        cont_cols = seq_along(contcols),
        lambda = 0,
        contkernel = contkernel,
        nomkernel = nomkernel,
        ordkernel = ordkernel
      )
      pyx_cont <- pxy_list_cont$py_x
      nearest_neighbours_ratios <- apply(pyx_cont, 2, function(x) max(x) / max(x[-which.max(x)]))
      if (any(!is.finite(nearest_neighbours_ratios))) return(Inf)
      mean(nearest_neighbours_ratios)
    }
  }
  result <- tryCatch({
    root <- uniroot(
      f = function(s) compute_avg(s) - ratio_thresh,
      interval = c(s_min, s_max),
      tol = 0.05
    )
    floor(root$root * 10) / 10
  }, error = function(e) {
    if (compute_avg(s_max) >= ratio_thresh) s_max else s_min
  })
  return(result)
}


find_s_fn <- function(X, contcols, catcols, contkernel, nomkernel, ordkernel,
                      cat_rel_imp, n_bw_reps = 100, s_min = 0.1, s_max = 100,
                      nystrom) {
  thresh <- (cat_rel_imp)^(length(catcols))
  if (nystrom) {
    n_bw_sample <- min(nrow(X), 1000)
    compute_avg <- function(s_val) {
      avg_py_x_reps <- replicate(n_bw_reps, {
        bw_sample_idx <- sample(nrow(X), n_bw_sample)
        X_bw <- X[bw_sample_idx, ]
        pxy_list_cont <- coord_to_pxy_R(
          as.data.frame(X_bw[, contcols]),
          s = s_val,
          cat_cols = c(),
          cont_cols = seq_along(contcols),
          lambda = 0,
          contkernel = contkernel,
          nomkernel = nomkernel,
          ordkernel = ordkernel
        )
        pyx_cont <- pxy_list_cont$py_x
        furthest_neighbours_ratios <- apply(pyx_cont, 2, function(x) max(x) / min(x[x > 0]))
        mean(furthest_neighbours_ratios)
      })
      if (all(!is.finite(avg_py_x_reps))) return(Inf)
      median(avg_py_x_reps[is.finite(avg_py_x_reps)])
    }
  } else {
    compute_avg <- function(s_val) {
      pxy_list_cont <- coord_to_pxy_R(
        as.data.frame(X[, contcols]),
        s = s_val,
        cat_cols = c(),
        cont_cols = seq_along(contcols),
        lambda = 0,
        contkernel = contkernel,
        nomkernel = nomkernel,
        ordkernel = ordkernel
      )
      pyx_cont <- pxy_list_cont$py_x
      avg_max_min_py_x <- mean(apply(pyx_cont, 2, function(x) max(x) / min(x[x > 0])))
      if (!is.finite(avg_max_min_py_x)) return(Inf)
      avg_max_min_py_x
    }
  }
  result <- tryCatch({
    root <- uniroot(
      f = function(s) compute_avg(s) - thresh,
      interval = c(s_min, s_max),
      tol = 0.05
    )
    s_candidate <- ceiling(root$root * 10) / 10
    while (s_candidate > s_min && compute_avg(s_candidate - 0.1) < thresh) {
      s_candidate <- s_candidate - 0.1
    }
    s_candidate
  }, error = function(e) {
    if (compute_avg(s_min) < thresh) s_min else s_max
  })
  return(result)
}

# Helper function to preprocess continuous data
preprocess_cont_data <- function(X) {
  X <- data.frame(X)
  X <- scale(X) 
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
compute_bandwidth_cont <- function(X, contcols, contkernel,
                                   nomkernel, ordkernel, nystrom){
  n_bw_reps <- 100
  s_val <- find_s_nn(X, contcols, contkernel, nomkernel, ordkernel,
                     n_bw_reps = n_bw_reps, s_min = 0.1, s_max = 100,
                     nystrom)
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
                             cat_first, nystrom){
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
    if (length(s) == 1 && s == -1){
      n_bw_reps <- 100
      s <- find_s_nn(X, contcols, contkernel, nomkernel, ordkernel,
                     n_bw_reps = n_bw_reps, s_min = 0.1, s_max = 100,
                     nystrom)
      n_bw_sample <- min(nrow(X), 1000)
      avg_py_x_reps <- replicate(n_bw_reps, {
        bw_sample_idx <- sample(nrow(X), n_bw_sample)
        X_bw <- X[bw_sample_idx, ]
        pxy_list_cont <- coord_to_pxy_R(
          as.data.frame(X_bw[, contcols]),
          s = s,
          cat_cols = c(),
          cont_cols = c(1:length(contcols)),
          lambda = 0,
          contkernel = contkernel,
          nomkernel = nomkernel,
          ordkernel = ordkernel)
        pyx_cont <- pxy_list_cont$py_x
        furthest_neighbours_ratios <- apply(pyx_cont, 2, FUN = function(x) max(x)/min(x[which(x > 0)]))
        mean(furthest_neighbours_ratios)
      })
      avg_max_min_py_x <- median(avg_py_x_reps[is.finite(avg_py_x_reps)])
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
    s <- find_s_fn(X, contcols, catcols, contkernel, nomkernel, ordkernel,
                   cat_rel_imp, n_bw_reps = 100, s_min = 0.1, s_max = 100,
                   nystrom)
  }
  bw <- rep(NA, ncol(X))
  bw[contcols] <- s
  bw[catcols] <- lambda
  bws_vec <- bw
  return(bws_vec)
}

