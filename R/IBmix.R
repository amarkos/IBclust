IBmix <- function(X, ncl, beta, randinit = NULL,
                  s = -1, lambda = -1, scale = TRUE,
                  maxiter = 100, nstart = 100, conv_tol = 1e-5,
                  contkernel = "gaussian",
                  nomkernel = "aitchisonaitken", ordkernel = "liracine",
                  cat_first = FALSE, verbose = FALSE) {
  
  # Validate inputs
  if (!is.numeric(ncl) || ncl <= 1 || ncl != round(ncl)) {
    stop("Input 'ncl' must be a positive integer greater than 1.")
  }
  if (!is.numeric(maxiter) || maxiter <= 0 || maxiter != round(maxiter)) {
    stop("'maxiter' must be a positive integer.")
  }
  if (!is.numeric(beta) || beta <= 0) {
    stop("Input 'beta' must be a positive number.")
  }
  if (!is.numeric(conv_tol) || conv_tol <= 0 || conv_tol >= 1) {
    stop("Input 'conv_tol' must be between 0 and 1.")
  }
  if (!is.numeric(nstart) || nstart <= 0 || nstart != round(nstart)) {
    stop("'nstart' must be a positive integer.")
  }
  if (!is.null(randinit) && (!is.numeric(randinit) || length(randinit) != nrow(X))) {
    stop("'randinit' must be a numeric vector with length equal to the number of rows in 'X', or NULL.")
  }
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
  
  py_x <- pxy_list$py_x
  px <- pxy_list$px
  pxy <- pxy_list$pxy
  hy <- pxy_list$hy
  
  ######################################################
  best_clust <- IBmix_iterate(X, ncl = ncl, beta = beta,
                              randinit = randinit, conv_tol,
                              tol = 0, py_x, hy, px, maxiter,
                              bws_vec, contcols, catcols,
                              runs = nstart, verbose = verbose)
  ######################################################
  
  # Wrap into an S3 object of class gibclust
  res <- new_gibclust(
    cluster = best_clust$Cluster,
    entropy = best_clust$Entropy,
    cond_entropy = best_clust$CondEntropy,
    mutual_info = best_clust$MutualInfo,
    info_xt = best_clust$InfoXT,
    beta = best_clust$beta,
    alpha = best_clust$alpha,
    s = best_clust$s,
    lambda = best_clust$lambda,
    call = match.call(),
    ncl = ncl,
    n = nrow(X),
    iters = ifelse(best_clust$converged,
                   as.integer(best_clust$iters),
                   maxiter),
    converged = best_clust$converged,
    conv_tol = conv_tol,
    contcols = contcols,
    catcols = catcols,
    kernels = list(cont = contkernel,
                   nom = nomkernel,
                   ord = ordkernel)
  )
  return(res)
}
