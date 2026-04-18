#' Information Bottleneck Clustering for Mixed-Type Data
#'
#' The \code{IBmix} function implements the Information Bottleneck (IB) algorithm
#' for clustering datasets containing continuous, categorical (nominal and ordinal), and mixed-type variables.
#' This method optimizes an information-theoretic objective to preserve
#' relevant information in the cluster assignments while achieving effective data compression
#' \insertCite{strouse_ib_2019}{IBclust}.
#'
#' @param X A data frame containing the input data to be clustered. It should include categorical variables
#'   (\code{factor} for nominal and \code{ordered} for ordinal) and continuous variables (\code{numeric}).
#' @param ncl An integer specifying the number of clusters.
#' @param beta Regularisation strength characterizing the tradeoff between compression and relevance. Must be non-negative.
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
#' @param conv_tol Convergence tolerance level; for a cluster membership matrix \eqn{U^{(m)}} at iteration \eqn{m}, convergence is achieved if
#'   \eqn{\sum_{i,j}\lvert U_{i,j}^{m+1} - U_{i,j}^m \rvert \le } \code{conv_tol}. Must be in range \eqn{[0, 1]}. Defaults to \code{1e-5}.
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
#'
#' @return An object of class \code{"gibclust"} representing the final clustering result. The returned object is a list with the following components:
#'   \item{Cluster}{An integer vector giving the cluster assignments for each data point.}
#'   \item{Entropy}{A numeric value representing the entropy of the cluster assignments at convergence.}
#'   \item{CondEntropy}{A numeric value representing the conditional entropy of cluster assignment, given the observation weights \eqn{H(T \mid X)}.}
#'   \item{MutualInfo}{A numeric value representing the mutual information, \eqn{I(Y;T)}, between the original labels (\eqn{Y}) and the cluster assignments (\eqn{T}).}
#'   \item{InfoXT}{A numeric value representing the mutual information, \eqn{I(X;T)}, between the original observations weights (\eqn{X}) and the cluster assignments (\eqn{T}).}
#'   \item{beta}{A numeric vector of the final beta values used in the iterative procedure.}
#'   \item{alpha}{A numeric value of the strength of conditional entropy used, controlling fuzziness of the solution. This is by default equal to \eqn{1} for \code{IBmix}.}
#'   \item{s}{A numeric vector of bandwidth parameters used for the continuous variables. A value of \eqn{-1} is returned if all variables are categorical.}
#'   \item{lambda}{A numeric vector of bandwidth parameters used for the categorical variables. A value of \eqn{-1} is returned if all variables are continuous.}
#'   \item{call}{The matched call.}
#'   \item{ncl}{Number of clusters.}
#'   \item{n}{Number of observations.}
#'   \item{iters}{Number of iterations used to obtain the returned solution.}
#'   \item{converged}{Logical indicating whether convergence was reached before \code{maxiter}.}
#'   \item{conv_tol}{Numeric convergence tolerance.}
#'   \item{contcols}{Indices of continuous columns in \code{X}.}
#'   \item{catcols}{Indices of categorical columns in \code{X}.}
#'   \item{kernels}{List with names of kernels used for continuous, nominal, and ordinal features.}
#'
#' Objects of class \code{"gibclust"} support the following methods:
#'   \itemize{
#'     \item \code{print.gibclust}: Display a concise description of the cluster assignment.
#'     \item \code{summary.gibclust}: Show detailed information including cluster sizes,
#'           information-theoretic metrics, hyperparameters, and convergence details.
#'     \item \code{plot.gibclust}: Produce diagnostic plots:
#'       \itemize{
#'         \item \code{type = "sizes"}: barplot of cluster sizes or hardened sizes (IB/GIB).
#'         \item \code{type = "info"}: barplot of entropy, conditional entropy, and mutual information.
#'         \item \code{type = "beta"}: trajectory of \eqn{\log \beta} over iterations (only available for hard clustering outputs obtained using \code{DIBmix}).
#'     }
#'   }
#'
#' @details
#' The \code{IBmix} function produces a fuzzy clustering of the data while retaining maximal information about the original variable
#' distributions. The Information Bottleneck algorithm optimizes an information-theoretic
#' objective that balances information preservation and compression. Bandwidth parameters for categorical
#' (nominal, ordinal) and continuous variables are adaptively determined if not provided. This iterative
#' process identifies stable and interpretable cluster assignments by maximizing mutual information while
#' controlling complexity. The method is well-suited for datasets with mixed-type variables and integrates
#' information from all variable types effectively. For data sets with over a thousand observations (\eqn{n > 1000}),
#' a Nystr\enc{ö}{o}m approximation of the kernel Gram matrix can be used for a quicker implementation \insertCite{williams2000using}{IBclust}.
#'
#' The following kernel functions can be used to estimate densities for the clustering procedure. For continuous variables:
#'
#' \itemize{
#'   \item \emph{Gaussian (RBF) kernel \insertCite{silverman_density_1998}{IBclust}:}
#'   \deqn{K_c\left(\frac{x - x'}{s}\right) = \frac{1}{\sqrt{2\pi}} \exp\left\{-\frac{\left(x - x'\right)^2}{2s^2}\right\}, \quad s > 0.}
#'   \item \emph{Epanechnikov kernel \insertCite{epanechnikov1969non}{IBclust}:}
#'   \deqn{K_c(x - x'; s) = \begin{cases}
#'     \frac{3}{4\sqrt{5}}\left(1 - \frac{(x-x')^2}{5s^2} \right), & \text{if } \frac{(x - x')^2}{s^2} < 5 \\
#'     0, & \text{otherwise}
#' \end{cases}, \quad s > 0.}
#' }
#'
#' For nominal (unordered categorical variables):
#'
#' \itemize{
#' \item \emph{Aitchison & Aitken kernel \insertCite{aitchison_kernel_1976}{IBclust}:}
#' \deqn{K_u(x = x'; \lambda) = \begin{cases}
#'     1 - \lambda, & \text{if } x = x' \\
#'     \frac{\lambda}{\ell - 1}, & \text{otherwise}
#' \end{cases}, \quad 0 \leq \lambda \leq \frac{\ell - 1}{\ell}.}
#' \item \emph{Li & Racine kernel \insertCite{ouyang2006cross}{IBclust}:}
#' \deqn{K_u(x = x'; \lambda) = \begin{cases}
#'     1, & \text{if } x = x' \\
#'     \lambda, & \text{otherwise}
#' \end{cases}, \quad 0 \leq \lambda \leq 1.}
#' }
#'
#' For ordinal (ordered categorical) variables:
#'
#' \itemize{
#' \item \emph{Li & Racine kernel \insertCite{li_nonparametric_2003}{IBclust}:}
#' \deqn{K_o(x = x'; \nu) = \begin{cases}
#'     1, & \text{if } x = x' \\
#'     \nu^{|x - x'|}, & \text{otherwise}
#' \end{cases}, \quad 0 \leq \nu \leq 1.}
#' \item \emph{Wang & van Ryzin kernel \insertCite{wang1981class}{IBclust}:}
#' \deqn{K_o(x = x'; \nu) = \begin{cases}
#'     1 - \nu, & \text{if } x = x' \\
#'     \frac{1-\nu}{2}\nu^{|x - x'|}, & \text{otherwise}
#' \end{cases}, \quad 0 \leq \nu \leq 1.}
#' }
#'
#' The bandwidth parameters \eqn{s}, \eqn{\lambda}, and \eqn{\nu} control the smoothness of the density estimate and are automatically determined by the algorithm if not provided by the user using the approach in \insertCite{costa_dib_2025;textual}{IBclust}. \eqn{\ell} is the number of levels of the categorical variable. For ordinal variables, the lambda parameter of the function is used to define \eqn{\nu}.
#'
#' @examples
#' # Example dataset with categorical, ordinal, and continuous variables
#' set.seed(123)
#' data_mix <- data.frame(
#'   cat_var = factor(sample(letters[1:3], 100, replace = TRUE)),      # Nominal categorical variable
#'   ord_var = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
#'                    levels = c("low", "medium", "high"),
#'                    ordered = TRUE),                                # Ordinal variable
#'   cont_var1 = rnorm(100),                                          # Continuous variable 1
#'   cont_var2 = runif(100)                                           # Continuous variable 2
#' )
#'
#' # Perform Mixed-Type Fuzzy Clustering
#' result_mix <- IBmix(X = data_mix, ncl = 3, beta = 2, nstart = 1)
#'
#' # Print clustering results
#' print(result_mix$Cluster)       # Cluster membership matrix
#' print(result_mix$InfoXT)        # Mutual information between X and T
#' print(result_mix$MutualInfo)    # Mutual information between Y and T
#'
#' # Summary of output
#' summary(result_mix)
#'
#' # Simulated categorical data example
#' set.seed(123)
#' data_cat <- data.frame(
#'   Var1 = as.factor(sample(letters[1:3], 100, replace = TRUE)),  # Nominal variable
#'   Var2 = as.factor(sample(letters[4:6], 100, replace = TRUE)),  # Nominal variable
#'   Var3 = factor(sample(c("low", "medium", "high"), 100, replace = TRUE),
#'                 levels = c("low", "medium", "high"), ordered = TRUE)  # Ordinal variable
#' )
#'
#' # Perform fuzzy clustering on categorical data with standard IB
#' result_cat <- IBmix(X = data_cat, ncl = 3, beta = 15, lambda = -1, nstart = 2, maxiter = 200)
#'
#' # Print clustering results
#' print(result_cat$Cluster)       # Cluster membership matrix
#' print(result_cat$InfoXT)        # Mutual information between X and T
#' print(result_cat$MutualInfo)    # Mutual information between Y and T
#'
#' plot(result_cat, type = "sizes") # Bar plot of cluster sizes (hardened assignments)
#' plot(result_cat, type = "info")  # Information-theoretic quantities plot
#'
#' # Simulated continuous data example
#' set.seed(123)
#' # Continuous data with 100 observations, 5 features
#' data_cont <- as.data.frame(matrix(rnorm(500), ncol = 5))
#'
#' # Perform fuzzy clustering on continuous data with standard IB
#' result_cont <- IBmix(X = data_cont, ncl = 3, beta = 50, s = -1, nstart = 2)
#'
#' # Print clustering results
#' print(result_cont$Cluster)       # Cluster membership matrix
#' print(result_cont$InfoXT)        # Mutual information between X and T
#' print(result_cont$MutualInfo)    # Mutual information between Y and T
#'
#' @author Efthymios Costa, Ioanna Papatsouma, Angelos Markos
#'
#' @references
#' \insertRef{strouse_ib_2019}{IBclust}
#'
#' \insertRef{aitchison_kernel_1976}{IBclust}
#'
#' \insertRef{li_nonparametric_2003}{IBclust}
#'
#' \insertRef{silverman_density_1998}{IBclust}
#'
#' \insertRef{ouyang2006cross}{IBclust}
#'
#' \insertRef{wang1981class}{IBclust}
#'
#' \insertRef{epanechnikov1969non}{IBclust}
#'
#' \insertRef{williams2000using}{IBclust}
#'
#' @keywords clustering
#' @export
#' @rdname IBmix
IBmix <- function(X, ncl, beta, randinit = NULL,
                  s = -1, lambda = -1, scale = TRUE,
                  maxiter = 100, nstart = 100, conv_tol = 1e-5,
                  contkernel = "gaussian",
                  nomkernel = "aitchisonaitken", ordkernel = "liracine",
                  cat_first = FALSE, verbose = FALSE, nystrom = FALSE,
                  n_landmarks = NULL) {
  
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
  if (!is.logical(nystrom)) {
    stop("'nystrom' must be a logical (TRUE or FALSE).")
  }
  if (nystrom){
    if (is.null(n_landmarks)){
      n_landmarks <- ceiling(sqrt(nrow(X)))
    }
  }
  prep_list <- input_checks_preprocess(X, s, lambda,
                                       scale, contkernel, nomkernel,
                                       ordkernel, cat_first,
                                       nystrom = nystrom,
                                       n_landmarks = n_landmarks)
  X <- prep_list$X
  bws_vec <- prep_list$bws_vec
  contcols <- prep_list$contcols
  catcols <- prep_list$catcols
  
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
                                       n_landmarks = n_landmarks)
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
  }
  
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
