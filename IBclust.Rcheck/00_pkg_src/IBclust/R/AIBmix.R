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
  obs_names <- rownames(X)
  if (is.null(obs_names)) obs_names <- as.character(seq_len(nrow(X)))
  
  res <- new_aibclust(
    merges = best_clust$merges,
    merge_costs = best_clust$merge_costs,
    partitions = best_clust$partitions,
    I_T_Y = best_clust$I_Z_Y,
    I_X_Y = best_clust$I_X_Y,
    info_ret = best_clust$info_ret,
    s = if (length(contcols) == 0) -1 else as.vector(bws_vec[contcols]),
    lambda = if (length(catcols) == 0) -1 else as.vector(bws_vec[catcols]),
    call = match.call(),
    n = nrow(X),
    contcols = contcols,
    catcols = catcols,
    kernels = list(cont = contkernel, nom = nomkernel, ord = ordkernel),
    obs_names = obs_names
  )
  return(res)
}
