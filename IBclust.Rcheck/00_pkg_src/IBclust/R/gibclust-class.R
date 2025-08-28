# Internal constructor for gibclust class objects

#' @keywords internal
#' @noRd
new_gibclust <- function(
    cluster, entropy, cond_entropy, mutual_info, info_xt,
    beta, alpha, s, lambda,
    call, ncl, n, iters = NA_integer_, converged = NA,
    conv_tol = NA_real_, contcols = integer(), catcols = integer(),
    kernels = list(cont = NA_character_,
                   nom = NA_character_,
                   ord = NA_character_)
) {
  cl <- cluster  # use local variable 'cl' to avoid confusion
  
  if (is.matrix(cl)) {
    d <- dim(cl)
    # If it's n x ncl, transpose to ncl x n
    if (length(d) == 2 && d[1] == n && d[2] == ncl) {
      cl <- t(cl)
    }
    storage.mode(cl) <- "double"
  } else {
    cl <- as.integer(cl)
  }
  
  x <- list(
    Cluster = cl,
    Entropy = as.numeric(entropy),
    CondEntropy = as.numeric(cond_entropy),
    MutualInfo = as.numeric(mutual_info),
    InfoXT = info_xt,
    beta = beta,
    alpha = alpha,
    s = s,
    lambda = lambda,
    call = call,
    ncl = as.integer(ncl),
    n = as.integer(n),
    iters = as.integer(iters),
    converged = isTRUE(converged),
    conv_tol = conv_tol,
    contcols = contcols,
    catcols = catcols,
    kernels = kernels
  )
  validate_gibclust(x)
  class(x) <- "gibclust"
  x
}

#' @keywords internal
#' @noRd
validate_gibclust <- function(x) {
  stopifnot(
    is.list(x),
    !is.null(x$n), !is.null(x$ncl),
    x$n >= 1L, x$ncl >= 1L
  )
  
  is_hard <- is.atomic(x$Cluster) && length(x$Cluster) == x$n
  is_soft <- is.matrix(x$Cluster) &&
    nrow(x$Cluster) == x$ncl &&
    ncol(x$Cluster) == x$n
  
  if (!is_hard && !is_soft) {
     stop("Cluster must be length-n vector (DIB) or (ncl x n) matrix (IB/GIB).")
  }
  
  if (is_soft) {
    cs <- colSums(x$Cluster)
    if (!all(is.finite(cs))) {
      stop("Non-finite values in membership matrix.")
    }
    # Allow a small tolerance for numerical noise
    if (!all(abs(cs - 1) < 1e-6)) {
      stop("Each column of Cluster (IB/GIB) must sum to 1.")
    }
  }
  stopifnot(
    length(x$MutualInfo) == 1L, is.numeric(x$MutualInfo),
    length(x$Entropy) == 1L, is.numeric(x$Entropy),
    length(x$CondEntropy)== 1L, is.numeric(x$CondEntropy),
    length(x$converged) == 1L || is.na(x$converged)
  )
  x
}