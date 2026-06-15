#' Sequential Information Bottleneck Clustering for Mixed-Type Data
#'
#' The \code{sIBmix} function implements the sequential Information Bottleneck
#' (sIB) algorithm for clustering
#' datasets containing continuous, categorical (nominal and ordinal), and
#' mixed-type variables. Unlike the agglomerative method (\code{\link{AIBmix}}),
#' sIB maintains exactly \code{ncl} clusters throughout and is guaranteed to
#' converge to a local maximum of the mutual information \eqn{I(Y; T)} \insertCite{slonim2002unsupervised}{IBclust}.
#'
#' @inheritParams DIBmix
#' @param eps Convergence tolerance. A sweep producing at most \code{eps * n}
#'   re-assignments declares convergence. Defaults to \eqn{0} (a fully stable
#'   sweep), consistent with the hard-clustering convention of \code{DIBmix}.
#'
#' @return An object of class \code{sibclust}. See
#'   \code{\link{sibclust-methods}} for the available S3 methods
#'   (\code{print}, \code{summary}, \code{plot}, \code{fitted},
#'   \code{coef}, \code{info_metrics}, \code{predict}).
#'
#' @details
#' The \code{sIBmix} function clusters data by maximising the information that the cluster assignments
#' retain about the original variable distributions, while holding the number of clusters fixed at \code{ncl}.
#' In contrast to the hierarchical agglomerative approach of \code{\link{AIBmix}}, the algorithm maintains a
#' partition into exactly \code{ncl} clusters throughout. At each step a single observation is removed from its
#' current cluster and re-assigned to the cluster for which the resulting loss of mutual information is smallest.
#' Each such move can only increase the retained information, or leave it unchanged, so the procedure is guaranteed
#' to converge to a local optimum. The criterion used to score each re-assignment is the same information-loss
#' measure that underlies \code{\link{AIBmix}}, so the two methods optimise the same objective through different
#' search strategies. Bandwidth parameters for categorical (nominal, ordinal) and continuous variables are
#' adaptively determined if not provided. The algorithm is run from several random initial partitions (controlled
#' by \code{nstart}) and the most informative solution is returned. The method is well-suited for datasets with
#' mixed-type variables and integrates information from all variable types effectively.
#'
#' See \code{\link{IBclust-package}} for details on the available kernel
#' functions and their bandwidth parameters.
#'
#' @examples
#' set.seed(123)
#' data_mix <- data.frame(
#'   cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),
#'   ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
#'                    levels = c("low", "medium", "high"), ordered = TRUE),
#'   cont_var1 = rnorm(100),
#'   cont_var2 = runif(100)
#' )
#'
#' # Hard partitional clustering with sequential IB
#' result <- sIBmix(X = data_mix, ncl = 3, nstart = 5)
#'
#' # Print output and provide summary
#' print(result)
#' summary(result)
#' fitted(result)               # Hard cluster labels
#' coef(result)                 # Bandwidths
#' info_metrics(result)         # information-theoretic quantities
#'
#' plot(result, type = "sizes")      # Plot of cluster sizes
#' plot(result, type = "info")       # Plot of information-theoretic quantities
#' plot(result, type = "importance") # Variable importance plot
#' plot(result, type = "similarity") # Similarity matrix plot
#'
#' predict(result, newdata = data_mix[1:5, ])   # Predict integer labels for new data
#'
#' @author Efthymios Costa, Ioanna Papatsouma, Angelos Markos
#'
#' @references
#' \insertRef{slonim2002unsupervised}{IBclust}
#' 
#' @keywords clustering
#' @export
sIBmix <- function(X, ncl, randinit = NULL,
                   s = -1, lambda = -1, scale = TRUE,
                   maxiter = 100, nstart = 100, eps = 0,
                   contkernel = "gaussian", nomkernel = "aitchisonaitken",
                   ordkernel = "liracine", cat_first = FALSE, verbose = FALSE,
                   keep_data = TRUE) {
  
  if (!is.numeric(ncl) || ncl <= 1 || ncl != round(ncl)) {
    stop("Input 'ncl' must be a positive integer greater than 1.")
  }
  if (!is.numeric(maxiter) || maxiter <= 0 || maxiter != round(maxiter)) {
    stop("'maxiter' must be a positive integer.")
  }
  if (!is.numeric(nstart) || nstart <= 0 || nstart != round(nstart)) {
    stop("'nstart' must be a positive integer.")
  }
  if (!is.numeric(eps) || eps < 0) {
    stop("'eps' must be a non-negative number.")
  }
  if (!is.null(randinit) && (!is.numeric(randinit) || length(randinit) != nrow(X))) {
    stop("'randinit' must be a numeric vector with length nrow(X), or NULL.")
  }
  X_original <- X
  prep_list <- input_checks_preprocess(X, s, lambda,
                                       scale, contkernel, nomkernel,
                                       ordkernel, cat_first,
                                       nystrom = FALSE,
                                       n_landmarks = NULL,
                                       landmark_indices = NULL,
                                       nystrom_available = FALSE,
                                       keep_data = keep_data)
  X <- prep_list$X
  bws_vec <- prep_list$bws_vec
  contcols <- prep_list$contcols
  catcols <- prep_list$catcols
  
  pxy_list <- coord_to_pxy_R(as.data.frame(X),
                             s = if (length(contcols) > 0) bws_vec[contcols] else -1,
                             lambda = if (length(catcols)  > 0) bws_vec[catcols]  else -1,
                             cat_cols = catcols,
                             cont_cols = contcols,
                             contkernel = contkernel,
                             nomkernel = nomkernel,
                             ordkernel = ordkernel)
  py_x <- pxy_list$py_x
  px <- pxy_list$px
  hy <- pxy_list$hy
  best_clust <- sIBmix_iterate(X, ncl = ncl, randinit = randinit,
                               py_x = py_x, hy = hy, px = px,
                               maxiter = maxiter, eps = eps,
                               bws_vec = bws_vec, contcols = contcols,
                               catcols = catcols, runs = nstart,
                               verbose = verbose)
  res <- new_sibclust(
    cluster = best_clust$Cluster,
    entropy = best_clust$Entropy,
    cond_entropy = best_clust$CondEntropy,
    mutual_info = best_clust$MutualInfo,
    info_xt = best_clust$InfoXT,
    s = best_clust$s,
    lambda = best_clust$lambda,
    call = match.call(),
    ncl = ncl,
    n = nrow(X),
    iters = best_clust$iters,
    converged = best_clust$converged,
    eps = eps,
    maxiter = maxiter,
    contcols = contcols,
    catcols = catcols,
    kernels = list(cont = contkernel, nom = nomkernel, ord = ordkernel),
    scale = scale
  )
  if (isTRUE(keep_data)) res$training_data <- X_original
  return(res)
}