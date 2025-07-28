AIBmix <- function(X, catcols, contcols, lambda = -1,
                 s = -1, scale = TRUE,
                 contkernel = "gaussian",
                 nomkernel = "aitchisonaitken",
                 ordkernel = "liracine",
                 cat_first = FALSE) {
  
  # Validate inputs
  if (!is.data.frame(X)) {
    stop("Input 'X' must be a data frame.")
  }
  if (!all(catcols %in% seq_along(X))) {
    stop("Some 'catcols' indices are out of bounds or invalid.")
  }
  
  if (!all(contcols %in% seq_along(X))) {
    stop("Some 'contcols' indices are out of bounds or invalid.")
  }
  
  if (any(duplicated(c(catcols, contcols)))) {
    stop("'catcols' and 'contcols' must not overlap.")
  }
  
  if (!is.logical(scale)) {
    stop("'scale' must be a logical value (TRUE or FALSE).")
  }
  if (!is.logical(cat_first)) {
    stop("'cat_first' must be a logical value (TRUE or FALSE).")
  }
  # Check kernel types
  if (!contkernel %in% c("gaussian", "epanechnikov")){
    stop("'contkernel' can only be one of 'gaussian' or 'epanechnikov'")
  }
  if (!nomkernel %in% c("aitchisonaitken", "liracine")){
    stop("'nomkernel' can only be one of 'aitchisonaitken' or 'liracine'")
  }
  if (!ordkernel %in% c("liracine", "wangvanryzin")){
    stop("'ordkernel' can only be one of 'liracine' or 'wangvanryzin'")
  }
  
  if (cat_first & any(c(s, lambda) != -1)){
    stop("'cat_first' can only be TRUE when all bandwidths are determined by the algorithm (s = -1, lambda = -1).")
  }
  
  # Validate lambda
  if (!is.numeric(lambda) ||
      !(length(lambda) == 1 || length(lambda) == length(catcols)) ||
      any(lambda <= 0 & lambda != -1)) {
    stop("'lambda' must be either a single numeric value (-1 for automatic selection or a positive value) or a numeric vector with positive values matching the number of 'catcols'.")
  }
  
  # Additional check for maximum lambda value for nominal variables
  if (length(lambda) > 1 && length(lambda) == length(catcols)) {
    if (nomkernel == "liracine"){
      max_lambda <- sapply(catcols, function(col) {
        l <- length(unique(X[, col]))
        (l - 1) / l
      })
    } else {
      max_lambda <- 1
    }
    if (any(lambda > max_lambda)) {
      stop("'lambda' values for nominal variables must not exceed their maximum allowable value.")
    }
  }
  
  
  # Validate s
  if (!is.numeric(s) ||
      !(length(s) == 1 || length(s) == length(contcols)) ||
      any(s <= 0 & s != -1)) {
    stop("'s' must be either a single numeric value (-1 for automatic selection or a positive value) or a numeric vector with positive values matching the number of 'contcols'.")
  }
  
  X <- data.frame(X)
  X[, catcols] <- preprocess_cat_data(X[, catcols])
  if (scale){
    X[, contcols] <- as.data.frame(preprocess_cont_data(X[, contcols]))
  }
  
  bws_vec <- compute_s_lambda(X, contcols, catcols, s, lambda,
                              contkernel, nomkernel, ordkernel,
                              cat_first)
  
  # Construct joint density with final bandwidths
  pxy_list <- coord_to_pxy_R(X, s = bws_vec[contcols],
                             cat_cols = catcols, cont_cols = contcols,
                             lambda = bws_vec[catcols],
                             contkernel = contkernel,
                             nomkernel = nomkernel,
                             ordkernel = ordkernel)
  
  pxy <- pxy_list$pxy
  # For AIB, we need strictly non-zero values in pxy...
  while (sum(pxy == 0) > 0){
    bws_vec[contcols] <- bws_vec[contcols] + 1e-1
    pxy_list <- coord_to_pxy_R(X, s = bws_vec[contcols],
                               cat_cols = catcols, cont_cols = contcols,
                               lambda = bws_vec[catcols],
                               contkernel = contkernel,
                               nomkernel = nomkernel,
                               ordkernel = ordkernel)
    pxy <- pxy_list$pxy
  }
  
  # Run AIB for hierarchical clustering
  best_clust <- AIB(pxy)
  dendrogram <- make_dendrogram(best_clust$merges,
                                best_clust$merge_costs,
                                labels = row.names(X))
  attr(dendrogram, "call") <- NULL
  best_clust[[length(best_clust)+1]] <- dendrogram
  names(best_clust)[length(best_clust)] <- "dendrogram"
  return(best_clust)
}
