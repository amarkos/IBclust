#' @keywords internal
#' @noRd
new_aibclust <- function(
    merges,
    merge_costs,
    partitions,
    I_T_Y,
    I_X_Y,
    info_ret,
    s,
    lambda,
    call,
    n,
    contcols = integer(),
    catcols  = integer(),
    kernels  = list(cont = NA_character_,
                    nom = NA_character_,
                    ord = NA_character_),
    obs_names = NULL
) {
  if (is.null(obs_names)) {
    obs_names <- as.character(seq_len(n))
  }
  x <- list(
    merges = merges,
    merge_costs = merge_costs,
    partitions = partitions,
    I_T_Y = I_T_Y,
    I_X_Y = I_X_Y,
    info_ret = info_ret,
    s = s,
    lambda = lambda,
    call = call,
    n = as.integer(n),
    contcols = as.integer(contcols),
    catcols = as.integer(catcols),
    kernels = kernels,
    obs_names = as.character(obs_names)
  )
  validate_aibclust(x)
  class(x) <- "aibclust"
  x
}

#' @keywords internal
#' @noRd
validate_aibclust <- function(x) {
  stopifnot(is.list(x), length(x$n) == 1L, x$n >= 1L)
  n <- x$n
  
  # merges: (n-1) x 2
  stopifnot(is.matrix(x$merges), ncol(x$merges) == 2L, nrow(x$merges) == n - 1L)
  storage.mode(x$merges) <- "integer"
  
  # merge_costs: length n-1
  stopifnot(is.numeric(x$merge_costs), length(x$merge_costs) == n - 1L)
  
  # partitions: list length n; each length n
  stopifnot(is.list(x$partitions), length(x$partitions) == n)
  for (p in x$partitions) stopifnot(length(p) == n)
  
  # info vectors
  stopifnot(is.numeric(x$I_T_Y), length(x$I_T_Y) == n)
  stopifnot(is.numeric(x$I_X_Y), length(x$I_X_Y) == 1L)
  stopifnot(is.numeric(x$info_ret), length(x$info_ret) == n)
  
  # bandwidths
  stopifnot(is.numeric(x$s), is.numeric(x$lambda))
  stopifnot(is.integer(x$contcols), is.integer(x$catcols))
  stopifnot(is.list(x$kernels))
  for (nm in c("cont", "nom", "ord")) stopifnot(nm %in% names(x$kernels))
  
  # labels
  stopifnot(is.character(x$obs_names), length(x$obs_names) == n)
  
  x
}