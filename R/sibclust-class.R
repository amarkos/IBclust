# Internal constructor for sibclust class objects
#' @keywords internal
#' @noRd
new_sibclust <- function(
    cluster, entropy, cond_entropy, mutual_info, info_xt,
    s, lambda,
    call, ncl, n, iters = NA_integer_, converged = NA,
    eps = 0, maxiter = NA_integer_, contcols = integer(), catcols = integer(),
    kernels = list(cont = NA_character_,
                   nom = NA_character_,
                   ord = NA_character_),
    scale = TRUE
) {
  cl <- as.integer(cluster)
  
  x <- list(
    Cluster = cl,
    Entropy = as.numeric(entropy),
    CondEntropy = as.numeric(cond_entropy),
    MutualInfo = as.numeric(mutual_info),
    InfoXT = info_xt,
    s = s,
    lambda = lambda,
    call = call,
    ncl = as.integer(ncl),
    n = as.integer(n),
    iters = as.integer(iters),
    converged = isTRUE(converged),
    eps = eps,
    maxiter = as.integer(maxiter),
    contcols = contcols,
    catcols = catcols,
    kernels = kernels,
    scale = scale
  )
  validate_sibclust(x)
  class(x) <- "sibclust"
  x
}
#' @keywords internal
#' @noRd
validate_sibclust <- function(x) {
  stopifnot(
    is.list(x),
    !is.null(x$n), !is.null(x$ncl),
    x$n >= 1, x$ncl >= 1
  )
  # sIB produces a hard partition
  if (!(is.atomic(x$Cluster) && length(x$Cluster) == x$n)) {
    stop("Cluster must be a length-n vector (sIB produces a hard partition).")
  }
  stopifnot(
    length(x$MutualInfo) == 1, is.numeric(x$MutualInfo),
    length(x$Entropy) == 1, is.numeric(x$Entropy),
    length(x$CondEntropy) == 1, is.numeric(x$CondEntropy),
    length(x$converged) == 1 || is.na(x$converged)
  )
  x
}