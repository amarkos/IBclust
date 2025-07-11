AIBcat <- function(X, lambda = -1,
                   nomkernel = "aitchisonaitken", ordkernel = "liracine"){
  
  # Validate inputs
  if (!is.data.frame(X)) {
    stop("Input 'X' must be a data frame.")
  }
  if (!nomkernel %in% c("aitchisonaitken", "liracine")){
    stop("'nomkernel' can only be one of 'aitchisonaitken' or 'liracine'")
  }
  if (!ordkernel %in% c("liracine", "wangvanryzin")){
    stop("'ordkernel' can only be one of 'liracine' or 'wangvanryzin'")
  }
  
  # Validate lambda
  if (!is.numeric(lambda) ||
      !(length(lambda) == 1 || length(lambda) == ncol(X)) ||
      any(lambda <= 0 & lambda != -1)) {
    stop("'lambda' must be either a single numeric value (-1 for automatic selection or a positive value) or a numeric vector with positive values matching the number of variables.")
  }
  
  # Additional check for maximum lambda value for nominal variables
  if (length(lambda) > 1 && length(lambda) == ncol(X)) {
    if (nomkernel == "liracine"){
      max_lambda <- sapply(1:ncol(X), function(col) {
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
  
  
  # Preprocessing
  X <- preprocess_cat_data(X)
  
  if (length(lambda) == 1){
    if (lambda == -1){
      # Compute lambda for categorical data
      lambda <- compute_lambda_cat(X, nomkernel, ordkernel)
    }
  }
  
  bws_vec <- lambda
  
  # Compute joint probability density for categorical variables
  pxy_list <- coord_to_pxy_R(as.data.frame(X), s = 0, cat_cols = seq_len(ncol(X)),
                             cont_cols = c(), lambda = bws_vec,
                             nomkernel = nomkernel,
                             ordkernel = ordkernel)
  pxy <- pxy_list$pxy
  
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
