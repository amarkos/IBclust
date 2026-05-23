AIB <- function(p_xy){
  n_x <- nrow(p_xy)
  
  pb <- txtProgressBar(style = 3, min = 0, max = n_x - 1)
  on.exit(close(pb), add = TRUE)
  clusters <- lapply(1:n_x, function(i) list(
    indices = i,
    p_z = p_xy[i,]/sum(p_xy[i,]),
    prob = sum(p_xy[i,])
  ))
  active <- rep(TRUE, n_x)
  D <- make_IB_distmat(p_xy)
  
  merges <- matrix(NA_integer_, n_x-1, 2)
  merge_costs <- numeric(n_x-1)
  partitions <- vector("list", n_x)
  I_Z_Y <- numeric(n_x)
  H_T   <- numeric(n_x)
  H_T_X <- numeric(n_x)
  I_T_X <- numeric(n_x)
  
  p_y <- colSums(p_xy)
  px  <- rowSums(p_xy)
  
  partitions[[n_x]] <- seq_len(n_x)
  I_Z_Y[n_x] <- mutual_information(p_xy)
  H_T[n_x]   <- entropy(px)
  H_T_X[n_x] <- 0
  I_T_X[n_x] <- H_T[n_x]
  
  for(step in seq_len(n_x-1)) {
    act_ids <- which(active)
    subD <- D[act_ids, act_ids, drop=FALSE]
    diag(subD) <- Inf
    
    ij <- which(subD == min(subD), arr.ind=TRUE)[1, ]
    i_idx <- act_ids[ij[1]]
    j_idx <- act_ids[ij[2]]
    
    merges[step, ] <- c(i_idx, j_idx)
    merge_costs[step] <- D[i_idx, j_idx]
    
    # do the merge
    ci <- clusters[[i_idx]]
    cj <- clusters[[j_idx]]
    p_i <- ci$prob
    p_j <- cj$prob
    p_new <- p_i + p_j
    p_z_new <- (p_i * ci$p_z + p_j * cj$p_z) / p_new
    
    clusters[[i_idx]] <- list(
      indices = c(ci$indices, cj$indices),
      p_z = p_z_new,
      prob = p_new
    )
    active[j_idx] <- FALSE
    D[j_idx, ] <- D[, j_idx] <- Inf
    
    # update only row/col i_idx in D
    for(k in act_ids) if(k != i_idx){
      D[i_idx, k] <- D[k, i_idx] <-
        (clusters[[i_idx]]$prob + clusters[[k]]$prob) *
        js_divergence(clusters[[i_idx]]$p_z,
                      clusters[[k]]$p_z)
    }
    
    # record partition at this step
    part <- integer(n_x)
    cid <- 1
    for(idx in which(active)){
      part[ clusters[[idx]]$indices ] <- cid
      cid <- cid + 1
    }
    partitions[[n_x - step]] <- part
    
    # build p(Y|T) matrix and cluster marginals
    active_ids <- which(active)
    Pmat <- t(sapply(active_ids, function(i) clusters[[i]]$p_z))
    qt <- sapply(active_ids, function(i) clusters[[i]]$prob)
    
    # compute I(Z;Y)
    ratio <- log2(Pmat/matrix(p_y, nrow(Pmat), ncol(Pmat), byrow=TRUE))
    ratio[is.infinite(ratio)] <- 0
    I_ZY_cur <- sum(qt*rowSums(Pmat*ratio))
    I_Z_Y[n_x-step] <- I_ZY_cur
    
    # compute info metrics
    ncl_now <- length(active_ids)
    qt_x <- matrix(0, nrow = ncl_now, ncol = n_x)
    for (k in seq_len(ncl_now)) {
      qt_x[k, part == k] <- 1
    }
    ht  <- entropy(qt)
    ht_x <- as.numeric(crossprod(px, entropy(qt_x)))
    H_T[n_x - step] <- ht
    H_T_X[n_x - step] <- ht_x
    I_T_X[n_x - step] <- ht - ht_x
    
    setTxtProgressBar(pb, step)
  }
  
  I_X_Y <- I_Z_Y[n_x]
  # compute info retained
  info_ret <- I_Z_Y/I_X_Y
  
  res_list <- list(
    merges = merges,
    merge_costs = merge_costs,
    partitions = partitions,
    I_Z_Y = I_Z_Y,
    I_X_Y = I_X_Y,
    info_ret = info_ret,
    H_T = H_T,
    H_T_X = H_T_X,
    I_T_X = I_T_X
  )
  return(res_list)
}