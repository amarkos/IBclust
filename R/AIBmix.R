AIBmix <- function(X, s = -1, lambda = -1,
                   scale = TRUE,
                   contkernel = "gaussian",
                   nomkernel = "aitchisonaitken",
                   ordkernel = "liracine",
                   cat_first = FALSE) {
  
  prep_list <- input_checks_preprocess(X, s, lambda,
                                       scale, contkernel, nomkernel,
                                       ordkernel, cat_first)
  X <- prep_list$X
  bws_vec <- prep_list$bws
  contcols <- prep_list$contcols
  catcols <- prep_list$catcols
  
  # Construct joint density with final bandwidths
  pxy_list <- coord_to_pxy_R(as.data.frame(X),
                             s = if (length(contcols) > 0){
                               bws_vec[contcols]
                               } else {
                                 -1
                                 },
                             lambda = if (length(catcols) > 0){
                               bws_vec[catcols]
                             } else {
                               -1
                             },
                             cat_cols = catcols,
                             cont_cols = contcols,
                             contkernel = contkernel,
                             nomkernel = nomkernel,
                             ordkernel = ordkernel)
  
  pxy <- pxy_list$pxy
  # For AIB, we need strictly non-zero values in pxy...
  while (sum(pxy == 0) > 0){
    bws_vec[contcols] <- bws_vec[contcols] + 1e-1
    pxy_list <- coord_to_pxy_R(as.data.frame(X),
                               s = if (length(contcols) > 0){
                                 bws_vec[contcols]
                               } else {
                                 -1
                               },
                               lambda = if (length(catcols) > 0){
                                 bws_vec[catcols]
                               } else {
                                 -1
                               },
                               cat_cols = catcols,
                               cont_cols = contcols,
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
  best_clust[[length(best_clust) + 1]] <- if (length(contcols) == 0) -1 else as.vector(bws_vec[contcols])
  names(best_clust)[length(best_clust)] <- "s"
  best_clust[[length(best_clust) + 1]] <- if (length(catcols) == 0) -1 else as.vector(bws_vec[catcols])
  names(best_clust)[length(best_clust)] <- "lambda"
  return(best_clust)
}
