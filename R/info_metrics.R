#' Extract information-theoretic metrics from an IBclust fit
#'
#' Returns the entropy, conditional entropy, and mutual information quantities
#' computed by the chosen Information Bottleneck variant. Methods are provided
#' for both \code{gibclust} and \code{aibclust} objects, returning a parallel
#' set of quantities with class-specific extras.
#'
#' The shared quantities across both classes are:
#' \itemize{
#'   \item \code{H_T}: the entropy \eqn{H(T)} of the cluster assignment.
#'   \item \code{H_T_X}: the conditional entropy \eqn{H(T \mid X)} of the cluster
#'     assignment given the observation weights. This is \eqn{0} for hard
#'     partitions, since \eqn{T} is a deterministic function of \eqn{X}.
#'   \item \code{I_T_X}: the mutual information \eqn{I(T; X)} between the
#'     cluster assignment and the observation weights. Equals \code{H_T} for
#'     hard partitions.
#'   \item \code{I_T_Y}: the mutual information \eqn{I(T; Y)} between the
#'     cluster assignment and the joint distribution \eqn{Y}.
#' }
#'
#' For \code{gibclust} objects, each quantity is a single numeric value
#' corresponding to the fitted partition. For \code{aibclust} objects, each
#' quantity is a numeric vector indexed by the number of clusters \eqn{m},
#' from \eqn{m = 1} (a single cluster) to \eqn{m = n} (all singletons).
#' Supplying \code{ncl} extracts the scalar values at the chosen cut.
#'
#' The \code{aibclust} output additionally includes:
#' \itemize{
#'   \item \code{I_X_Y}: the scalar baseline \eqn{I(X; Y)}, the mutual
#'     information between observation weights and the joint distribution.
#'   \item \code{info_ret}: the proportion of baseline information retained,
#'     \eqn{I(T_m; Y) / I(X; Y)}. Either a vector over \eqn{m} or a scalar at
#'     the requested cut.
#' }
#'
#' @param object A \code{gibclust} or \code{aibclust} object.
#' @param ncl For \code{aibclust} objects, the number of clusters at which to
#'   evaluate the metrics. If \code{NULL} (default), returns the full vectors
#'   of metrics over all cluster counts. Ignored for \code{gibclust} objects.
#' @param ... Additional arguments.
#'
#' @return A named list of information-theoretic quantities. For
#'   \code{gibclust}, contains \code{H_T}, \code{H_T_X}, \code{I_T_X},
#'   \code{I_T_Y}. For \code{aibclust}, contains those four quantities plus
#'   \code{I_X_Y} and \code{info_ret}.
#'
#' @seealso \code{\link{DIBmix}}, \code{\link{IBmix}}, \code{\link{GIBmix}},
#'   \code{\link{AIBmix}}.
#'
#' @examples
#' # gibclust: single values per quantity
#' fit_dib <- DIBmix(iris[, -5], ncl = 3, nstart = 5)
#' info_metrics(fit_dib)
#'
#' # aibclust: full vectors over cluster counts
#' fit_aib <- AIBmix(iris[, -5])
#' metrics_all <- info_metrics(fit_aib)
#' str(metrics_all)
#'
#' # aibclust: scalar values at a specific cut
#' info_metrics(fit_aib, ncl = 3)
#'
#' @export
info_metrics <- function(object, ...) UseMethod("info_metrics")

#' @rdname info_metrics
#' @method info_metrics gibclust
#' @exportS3Method
info_metrics.gibclust <- function(object, ...) {
  list(
    H_T = object$Entropy,
    H_T_X = object$CondEntropy,
    I_T_X = object$InfoXT,
    I_T_Y = object$MutualInfo
  )
}

#' @rdname info_metrics
#' @method info_metrics aibclust
#' @exportS3Method
info_metrics.aibclust <- function(object, ncl = NULL, ...) {
  if (is.null(ncl)) {
    return(list(
      H_T = object$H_T,
      H_T_X = object$H_T_X,
      I_T_X = object$I_T_X,
      I_T_Y = object$I_T_Y,
      I_X_Y = object$I_X_Y,
      info_ret = object$info_ret
    ))
  }
  if (!is.numeric(ncl) || length(ncl) != 1 || ncl != round(ncl)) {
    stop("Number of clusters 'ncl' must be a single integer.")
  }
  ncl <- as.integer(ncl)
  if (ncl < 1 || ncl > length(object$I_T_Y)) {
    stop(sprintf("Number of clusters 'ncl' must be between 1 and %d.", length(object$I_T_Y)))
  }
  list(
    H_T = object$H_T[ncl],
    H_T_X = object$H_T_X[ncl],
    I_T_X = object$I_T_X[ncl],
    I_T_Y = object$I_T_Y[ncl],
    I_X_Y = object$I_X_Y,
    info_ret = object$info_ret[ncl]
  )
}