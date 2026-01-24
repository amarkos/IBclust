AIB_nystrom <- function(B, col_sums) {
  n_x <- nrow(B)
  
  # Precompute conditional distributions py_x for each observation
  p_y_x_list <- vector("list", n_x)
  row_probs <- numeric(n_x)
  
  for (i in seq_len(n_x)) {
    unnorm <- as.vector(B %*% B[i, ])
    if (col_sums[i] > 1e-10) {
      unnorm <- unnorm / col_sums[i]
    }
    total <- sum(unnorm)
    if (total > 1e-10) {
      p_y_x_list[[i]] <- unnorm / total
      row_probs[i] <- total
    } else {
      p_y_x_list[[i]] <- rep(0, n_x)
      p_y_x_list[[i]][i] <- 1
      row_probs[i] <- 1 / n_x
    }
  }
  row_probs <- row_probs / sum(row_probs)
  
  # Initial clusters
  clusters <- lapply(seq_len(n_x), function(i) list(
    indices = i,
    p_z = p_y_x_list[[i]],
    prob = row_probs[i]
  ))
  active <- rep(TRUE, n_x)
  
  # Precompute distance matrix
  D <- make_IB_distmat_nystrom_cpp(B, col_sums)
  
  merges <- matrix(NA_integer_, n_x - 1, 2)
  merge_costs <- numeric(n_x - 1)
  partitions <- vector("list", n_x)
  I_Z_Y <- numeric(n_x)
  
  # Record initial state
  partitions[[n_x]] <- seq_len(n_x)
  I_Z_Y[n_x] <- mutual_information_nystrom_cpp(B, col_sums)
  
  # Precompute p(y)
  p_y <- Reduce(`+`, lapply(seq_len(n_x), function(i) row_probs[i] * p_y_x_list[[i]]))
  
  for (step in seq_len(n_x - 1)) {
    act_ids <- which(active)
    subD <- D[act_ids, act_ids, drop = FALSE]
    diag(subD) <- Inf
    
    ij <- which(subD == min(subD), arr.ind = TRUE)[1, ]
    i_idx <- act_ids[ij[1]]
    j_idx <- act_ids[ij[2]]
    
    merges[step, ] <- c(i_idx, j_idx)
    merge_costs[step] <- D[i_idx, j_idx]
    
    # Merge
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
    
    # Update distances
    for (k in act_ids) {
      if (k != i_idx) {
        D[i_idx, k] <- D[k, i_idx] <-
          (clusters[[i_idx]]$prob + clusters[[k]]$prob) *
          js_divergence(clusters[[i_idx]]$p_z, clusters[[k]]$p_z)
      }
    }
    
    # Record partition
    part <- integer(n_x)
    cid <- 1
    for (idx in which(active)) {
      part[clusters[[idx]]$indices] <- cid
      cid <- cid + 1
    }
    partitions[[n_x - step]] <- part
    
    # Compute I(Z;Y)
    active_ids <- which(active)
    Pmat <- t(sapply(active_ids, function(i) clusters[[i]]$p_z))
    w <- sapply(active_ids, function(i) clusters[[i]]$prob)
    
    ratio <- log2(Pmat / matrix(p_y, nrow(Pmat), ncol(Pmat), byrow = TRUE))
    ratio[is.infinite(ratio)] <- 0
    
    I_ZY_cur <- sum(w * rowSums(Pmat * ratio))
    I_Z_Y[n_x - step] <- I_ZY_cur
  }
  
  I_X_Y <- I_Z_Y[n_x]
  info_ret <- I_Z_Y / I_X_Y
  
  res_list <- list(
    merges = merges,
    merge_costs = merge_costs,
    partitions = partitions,
    I_Z_Y = I_Z_Y,
    I_X_Y = I_X_Y,
    info_ret = info_ret
  )
  return(res_list)
}