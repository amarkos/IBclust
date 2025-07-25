AIBcont <- function(X, s = -1, scale = TRUE, contkernel = "gaussian"){
  # Validate inputs
  if (!is.logical(scale)) {
    stop("'scale' must be a logical value (TRUE or FALSE).")
  }
  # Check kernel types
  if (!contkernel %in% c("gaussian", "epanechnikov")){
    stop("'contkernel' can only be one of 'gaussian' or 'epanechnikov'")
  }
  # Validate s
  if (!is.numeric(s) ||
      !(length(s) == 1 || length(s) == ncol(X)) ||
      any(s <= 0 & s != -1)) {
    stop("'s' must be either a single numeric value (-1 for automatic selection or a positive value) or a numeric vector with positive values matching the number of 'contcols'.")
  }
  # Preprocessing
  if (scale){
    X <- as.data.frame(preprocess_cont_data(X))
  }
  
  # Bandwidth computation
  if (length(s) == 1){
    if (s == -1){
      s <- compute_bandwidth_cont(X, contkernel = contkernel)
    }
  }
  
  # Compute joint probability density for continuous variables
  pxy_list <- coord_to_pxy_R(as.data.frame(X), s = s, cat_cols = c(),
                             cont_cols = seq_len(ncol(X)), lambda = 0,
                             contkernel = contkernel)
  pxy <- pxy_list$pxy
  # For AIB, we need strictly non-zero values in pxy...
  while (sum(pxy == 0) > 0){
    s <- s + 1e-1
    pxy_list <- coord_to_pxy_R(as.data.frame(X), s = s, cat_cols = c(),
                               cont_cols = seq_len(ncol(X)), lambda = 0,
                               contkernel = contkernel)
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
