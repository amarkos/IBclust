#' Internal Function: sIBmix_iterate
#'
#' Performs sequential ("draw and re-assign") clustering for the sIBmix algorithm.
#'
#' @param X A data frame or matrix containing the dataset.
#' @param ncl Number of clusters.
#' @param randinit Optional initial cluster assignments.
#' @param eps Convergence tolerance; a sweep with at most \code{eps * n} re-assignments declares convergence.
#' @param py_x Conditional probability matrix \( p(y|x) \).
#' @param hy Entropy \( H(Y) \).
#' @param px Probability matrix \( p(x) \).
#' @param maxiter Maximum number of sweeps over the data.
#' @param bws_vec Bandwidth vector.
#' @param contcols Indices of continuous columns.
#' @param catcols Indices of categorical columns.
#' @param runs Number of random starts.
#' @param verbose Defaults to FALSE to suppress progress messages. Change to TRUE to print.
#'
#' @return A list containing clustering results.
#'
#' @keywords internal
#' @noRd
sIBmix_iterate <- function(X, ncl, randinit, eps, py_x, hy, px, maxiter,
                           bws_vec, contcols, catcols, runs, verbose = FALSE) {
  
  n  <- nrow(X)
  nY <- nrow(py_x)
  best_clust <- list(
    Cluster = rep(NA_integer_, n),
    Entropy = Inf,
    CondEntropy = Inf,
    MutualInfo = -Inf,
    InfoXT = Inf,
    s = if (length(contcols) == 0) -1 else as.vector(bws_vec[contcols]),
    lambda = if (length(catcols)  == 0) -1 else as.vector(bws_vec[catcols]),
    iters = NA_integer_,
    converged = NA
  )
  if (ncl == 1) {
    best_clust$Cluster <- rep(1, n)
    best_clust$Entropy <- 0
    best_clust$CondEntropy <- 0
    best_clust$MutualInfo <- 0
    best_clust$InfoXT <- 0
    best_clust$iters <- 0
    best_clust$converged <- FALSE
    return(best_clust)
  }
  
  H <- function(p) entropySingle(pmax(as.numeric(p), 0))
  
  pb <- txtProgressBar(style = 3, min = 0, max = runs)
  on.exit(close(pb), add = TRUE)
  
  for (run in seq_len(runs)) {
    setTxtProgressBar(pb, run)
    # Initial partition into ncl non-empty clusters
    if (is.null(randinit)) {
      assign <- sample(rep(1:ncl, each = ceiling(n / ncl)), size = n)
    } else {
      assign <- randinit
    }
    
    S <- matrix(0, nrow = nY, ncol = ncl)
    counts <- integer(ncl)
    for (k in seq_len(ncl)) {
      idx_k <- which(assign == k)
      counts[k] <- length(idx_k)
      S[, k] <- rowSums(py_x[, idx_k, drop = FALSE])
    }
    Hc <- vapply(seq_len(ncl), function(k) H(S[, k] / counts[k]), numeric(1))
    
    # Sweeps
    converged_run <- FALSE
    sweep_idx <- 0
    while (sweep_idx < maxiter) {
      sweep_idx <- sweep_idx + 1
      changes <- 0
      order_i <- sample.int(n)
      
      for (i in order_i) {
        a <- assign[i]
        if (counts[a] == 1) next
        
        col <- py_x[, i]
        
        # draw x out of cluster a
        S[, a] <- S[, a] - col
        counts[a] <- counts[a] - 1L
        Hc[a] <- H(S[, a] / counts[a])
        
        # exact information loss of re-inserting x into each cluster t
        score <- numeric(ncl)
        for (t in seq_len(ncl)) {
          nt <- counts[t]
          m_t <- (col + S[, t]) / (nt + 1)
          score[t] <- (nt + 1) * H(m_t) - nt * Hc[t]
        }
        b <- which.min(score)
        
        # assign x
        S[, b] <- S[, b] + col
        counts[b] <- counts[b] + 1
        Hc[b] <- H(S[, b] / counts[b])
        assign[i] <- b
        if (b != a) changes <- changes + 1
      }
      
      if (changes <= eps * n) {
        converged_run <- TRUE
        break
      }
    }
    
    # Metrics
    qt_x <- matrix(0, nrow = ncl, ncol = n)
    qt_x[cbind(assign, seq_len(n))] <- 1
    qt <- counts / n
    qy_t <- qy_t_step(py_x, qt_x, qt, px)
    metrics <- calc_metrics(beta = 1, qt, qy_t, hy, px, qt_x, quiet = TRUE)
    
    if (as.numeric(metrics$iyt) > best_clust$MutualInfo) {
      best_clust$Cluster <- assign
      best_clust$Entropy <- as.numeric(metrics$ht)
      best_clust$CondEntropy <- as.numeric(metrics$ht_x)
      best_clust$MutualInfo <- as.numeric(metrics$iyt)
      best_clust$InfoXT <- as.numeric(metrics$ixt)
      best_clust$iters <- sweep_idx
      best_clust$converged <- converged_run
    }
    if (verbose) message("Run ", run, " complete.\n")
  }
  
  return(best_clust)
}