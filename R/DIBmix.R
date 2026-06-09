#' Deterministic Information Bottleneck Clustering for Mixed-Type Data
#'
#' The \code{DIBmix} function implements the Deterministic Information Bottleneck (DIB) algorithm
#' for clustering datasets containing continuous, categorical (nominal and ordinal), and mixed-type variables.
#' This method optimises an information-theoretic objective to preserve
#' relevant information in the cluster assignments while achieving effective data compression
#' \insertCite{costa_dib_2025}{IBclust}.
#'
#' @param X A data frame containing the input data to be clustered. It should include categorical variables
#'   (\code{factor} for nominal and \code{ordered} for ordinal) and continuous variables (\code{numeric}).
#' @param ncl An integer specifying the number of clusters.
#' @param randinit An optional vector specifying the initial cluster assignments. If \code{NULL}, cluster assignments are initialized randomly.
#' @param s A numeric value or vector specifying the bandwidth parameter(s) for continuous variables. The values must be greater than \eqn{0}.
#'   The default value is \eqn{-1}, which enables the automatic selection of optimal bandwidth(s). Argument is ignored when no variables are continuous.
#' @param lambda A numeric value or vector specifying the bandwidth parameter for categorical variables. The default value is \eqn{-1}, which enables
#'   automatic determination of the optimal bandwidth. For nominal variables and \code{nomkernel = 'aitchisonaitken'}, the maximum allowable value of
#'   \code{lambda} is \eqn{(l - 1)/l}, where \eqn{l} represents the number of categories, whereas for \code{nomkernel = 'liracine'} the maximum
#'   allowable value is \eqn{1}. For ordinal variables, the maximum allowable value of \code{lambda} is \eqn{1}, regardless of what \code{ordkernel}
#'   is being used. Argument is ignored when all variables are continuous.
#' @param scale A logical value indicating whether the continuous variables should be scaled to have unit variance before clustering. Defaults to
#'   \code{TRUE}. Argument is ignored when all variables are categorical.
#' @param maxiter The maximum number of iterations allowed for the clustering algorithm. Defaults to \eqn{100}.
#' @param nstart The number of random initializations to run. The best clustering solution is returned. Defaults to \eqn{100}.
#' @param contkernel Kernel used for continuous variables. Can be one of \code{gaussian} (default) or \code{epanechnikov}. Argument is ignored when no
#'   variables are continuous.
#' @param nomkernel Kernel used for nominal (unordered categorical) variables. Can be one of \code{aitchisonaitken} (default) or \code{liracine}.
#'   Argument is ignored when no variables are nominal.
#' @param ordkernel Kernel used for ordinal (ordered categorical) variables. Can be one of \code{liracine} (default) or \code{wangvanryzin}. Argument
#'   is ignored when no variables are ordinal.
#' @param cat_first A logical value indicating whether bandwidth selection is prioritised for the categorical variables, instead of the continuous.
#'   Defaults to \code{FALSE}. Set to \code{TRUE} if you suspect that the continuous variables are not informative of the cluster structure. Can only
#'   be \code{TRUE} when data is of mixed-type and all bandwidths are selected automatically (i.e. \code{s = -1}, \code{lambda = -1}).
#' @param verbose Logical. Defaults to \code{FALSE} to suppress progress messages. Change to \code{TRUE} to print.
#' @param nystrom Logical. Indicates if the Nystr\enc{ö}{o}m approximation for kernel Gram matrices is to be used for quicker implementation. Defaults
#'   to \code{FALSE}. Change to \code{TRUE} only for data sets with over 1000 observations.
#' @param n_landmarks Number of randomly drawn landmark points used for the Nystr\enc{ö}{o}m approximation. Must be a positive integer less than the
#'   number of observations \code{nrow(X)}. Defaults to \code{NULL}, which selects \eqn{\lceil \sqrt{n} \rceil} observations, where \eqn{n} is the
#'   number of observations. Argument is ignored if \code{nystrom = FALSE}.
#' @param landmark_indices Optional integer vector specifying the exact indices
#'   of observations to use as landmark points for the Nystr\enc{ö}{o}m
#'   approximation. Must contain unique integers in \eqn{[1, n]}, where
#'   \eqn{n} is the number of observations. When provided, this overrides
#'   random landmark sampling; if \code{n_landmarks} is also supplied, its
#'   value must equal \code{length(landmark_indices)}. Defaults to \code{NULL},
#'   in which case landmarks are sampled randomly. Argument is ignored if
#'   \code{nystrom = FALSE}.
#' @param keep_data Logical; if \code{TRUE}, the original input data \code{X}
#'   is stored in the returned object as \code{training_data}, enabling the
#'   use of \code{predict()} and certain plotting methods without re-passing
#'   the data. Defaults to \code{FALSE} to keep returned objects lightweight.
#'
#' @return An object of class \code{"gibclust"}. See
#'   \code{\link{gibclust-methods}} for the available S3 methods
#'   (\code{print}, \code{summary}, \code{plot}, \code{fitted},
#'   \code{coef}, \code{info_metrics}, \code{predict}).
#'
#' @details
#' The \code{DIBmix} function clusters data while retaining maximal information about the original variable
#' distributions. The Deterministic Information Bottleneck algorithm optimises an information-theoretic
#' objective that balances information preservation and compression. Bandwidth parameters for categorical
#' (nominal, ordinal) and continuous variables are adaptively determined if not provided. This iterative
#' process identifies stable and interpretable cluster assignments by maximising mutual information while
#' controlling complexity. The method is well-suited for datasets with mixed-type variables and integrates
#' information from all variable types effectively. For data sets with over a thousand observations (\eqn{n > 1000}),
#' a Nystr\enc{ö}{o}m approximation of the kernel Gram matrix can be enabled via \code{nystrom = TRUE}; see \code{\link{IBclust-package}} for details.
#'
#' See \code{\link{IBclust-package}} for details on the available kernel
#' functions and their bandwidth parameters.
#'
#' @examples
#' # Example 1: Basic Mixed-Type Clustering
#' set.seed(123)
#'
#' # Create a more realistic dataset with mixed variable types
#' data_mix <- data.frame(
#'   # Categorical variables
#'   education = factor(sample(c("High School", "Bachelor", "Master", "PhD"), 150,
#'                            replace = TRUE, prob = c(0.4, 0.3, 0.2, 0.1))),
#'   employment = factor(sample(c("Full-time", "Part-time", "Unemployed"), 150,
#'                             replace = TRUE, prob = c(0.6, 0.25, 0.15))),
#'
#'   # Ordinal variable
#'   satisfaction = factor(sample(c("Low", "Medium", "High"), 150, replace = TRUE),
#'                        levels = c("Low", "Medium", "High"), ordered = TRUE),
#'
#'   # Continuous variables
#'   income = rlnorm(150, meanlog = 10, sdlog = 0.5),  # Log-normal income
#'   age = rnorm(150, mean = 35, sd = 10),             # Normally distributed age
#'   experience = rpois(150, lambda = 8)               # Years of experience
#' )
#'
#' # Perform Mixed-Type Clustering
#' result_mix <- DIBmix(X = data_mix, ncl = 3, nstart = 5)
#'
#' # View results
#' print(result_mix)
#' summary(result_mix)
#' table(fitted(result_mix))
#'
#' # Example 2: Comparing cat_first parameter
#' # When categorical variables are more informative
#' result_cat_first <- DIBmix(X = data_mix, ncl = 3,
#'                            cat_first = TRUE,  # Prioritise categorical variables
#'                            nstart = 5)
#'
#' # When continuous variables are more informative (default)
#' result_cont_first <- DIBmix(X = data_mix, ncl = 3,
#'                             cat_first = FALSE,
#'                             nstart = 5)
#'
#' # Compare clustering performance
#' if (requireNamespace("mclust", quietly = TRUE)){  # For adjustedRandIndex
#'   print(paste("Agreement between approaches:",
#'               round(mclust::adjustedRandIndex(fitted(result_cat_first),
#'                     fitted(result_cont_first)), 3)))
#'   }
#'
#' plot(result_cat_first, type = "sizes") # Bar plot of cluster sizes
#' plot(result_cat_first, type = "info")  # Information-theoretic quantities plot
#' plot(result_cat_first, type = "beta")  # Plot of evolution of beta
#' plot(result_cat_first, type = "importance") # Variable importance plot
#' plot(result_cat_first, type = "similarity") # Similarity matrix plot
#'
#' # Simulated categorical data example
#' data_cat <- data.frame(
#'   Var1 = as.factor(sample(letters[1:3], 200, replace = TRUE)),  # Nominal variable
#'   Var2 = as.factor(sample(letters[4:6], 200, replace = TRUE)),  # Nominal variable
#'   Var3 = factor(sample(c("low", "medium", "high"), 200, replace = TRUE),
#'                 levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
#' )
#'
#' # Perform hard clustering on categorical data with Deterministic IB
#' result_cat <- DIBmix(X = data_cat, ncl = 3, lambda = -1, nstart = 5)
#'
#' # Print clustering results
#' fitted(result_cat)        # Cluster assignments
#' info_metrics(result_cat)  # Information-theoretic quantities
#' coef(result_cat)          # Hyperparameter values
#'
#' # Simulated continuous data example
#' set.seed(123)
#' # Continuous data with 200 observations, 5 features
#' data_cont <- as.data.frame(matrix(rnorm(1000), ncol = 5))
#'
#' # Perform hard clustering on continuous data with Deterministic IB
#' result_cont <- DIBmix(X = data_cont, ncl = 3, s = -1, nstart = 5)
#'
#' # Print clustering results
#' fitted(result_cont)       # Cluster assignments
#'
#' # Summary of output
#' print(result_cont)
#' summary(result_cont)
#'
#' @author Efthymios Costa, Ioanna Papatsouma, Angelos Markos
#'
#' @references
#' \insertRef{costa_dib_2025}{IBclust}
#'
#' @keywords clustering
#' @export
#' @rdname DIBmix
DIBmix <- function(X, ncl, randinit = NULL,
                   s = -1, lambda = -1, scale = TRUE,
                   maxiter = 100, nstart = 100, contkernel = "gaussian",
                   nomkernel = "aitchisonaitken", ordkernel = "liracine",
                   cat_first = FALSE, verbose = FALSE, nystrom = FALSE,
                   n_landmarks = NULL, landmark_indices = NULL,
                   keep_data = TRUE) {

  # Validate inputs
  if (!is.numeric(ncl) || ncl <= 1 || ncl != round(ncl)) {
    stop("Input 'ncl' must be a positive integer greater than 1.")
  }
  if (!is.numeric(maxiter) || maxiter <= 0 || maxiter != round(maxiter)) {
    stop("'maxiter' must be a positive integer.")
  }
  if (!is.numeric(nstart) || nstart <= 0 || nstart != round(nstart)) {
    stop("'nstart' must be a positive integer.")
  }
  if (!is.null(randinit) && (!is.numeric(randinit) || length(randinit) != nrow(X))) {
    stop("'randinit' must be a numeric vector with length equal to the number of rows in 'X', or NULL.")
  }
  if (!is.logical(nystrom)) {
    stop("'nystrom' must be a logical (TRUE or FALSE).")
  }
  X_original <- X
  prep_list <- input_checks_preprocess(X, s, lambda,
                                       scale, contkernel, nomkernel,
                                       ordkernel, cat_first,
                                       nystrom = nystrom,
                                       n_landmarks = n_landmarks,
                                       landmark_indices = landmark_indices,
                                       keep_data = keep_data)
  X <- prep_list$X
  bws_vec <- prep_list$bws_vec
  contcols <- prep_list$contcols
  catcols <- prep_list$catcols
  n_landmarks <- prep_list$n_landmarks
  landmark_indices <- prep_list$landmark_indices
  
  # Construct joint density with final bandwidths
  if (nystrom){
    pxy_list <- coord_to_pxy_nystrom_R(as.data.frame(X),
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
                                       ordkernel = ordkernel,
                                       n_landmarks = n_landmarks,
                                       landmark_indices = landmark_indices)
    nystrom_landmarks <- pxy_list$landmark_indices
  } else {
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
    nystrom_landmarks <- NULL
  }
  py_x <- pxy_list$py_x
  px <- pxy_list$px
  pxy <- pxy_list$pxy
  hy <- pxy_list$hy

  ######################################################
  best_clust <- DIBmix_iterate(X, ncl = ncl, randinit = randinit,
                               tol = 0, py_x, hy, px, maxiter,
                               bws_vec, contcols, catcols,
                               runs = nstart, verbose = verbose)
  ######################################################
  if (best_clust$MutualInfo == Inf){
    warning("Initial cluster assignment remained unchanged; use other hyperparameter values for DIBmix to converge.")
  }
  
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
    conv_tol = 0,
    contcols = contcols,
    catcols = catcols,
    kernels = list(cont = contkernel,
                   nom = nomkernel,
                   ord = ordkernel),
    nystrom_landmarks = nystrom_landmarks,
    scale = scale
  )
  if (isTRUE(keep_data)) {
    res$training_data <- X_original
  }
  return(res)
}
