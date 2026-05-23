#' Extract information-theoretic metrics from an IBclust fit
#'
#' Returns the entropy, conditional entropy, and mutual information quantities
#' computed by the chosen IB variant. Methods are provided for both
#' \code{gibclust} and \code{aibclust} objects.
#'
#' For \code{gibclust} objects, a single value of each quantity is returned
#' corresponding to the fitted partition. For \code{aibclust} objects, the
#' quantities are returned as vectors indexed by the number of clusters
#' \eqn{m}, unless \code{ncl} is supplied (in which case scalar values at
#' the requested cluster count are returned). The \code{aibclust} output
#' additionally includes \code{I_X_Y} (a scalar baseline) and \code{info_ret}
#' (the fraction of baseline information retained at each cut).
#'
#' @param object A \code{gibclust} or \code{aibclust} object.
#' @param ncl For \code{aibclust} objects, the number of clusters at which to
#'   evaluate the metrics. If \code{NULL} (default), returns the full vector
#'   of metrics over all cluster counts. Ignored for \code{gibclust} objects.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A named list of information-theoretic quantities.
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
  if (!is.numeric(ncl) || length(ncl) != 1L || ncl != round(ncl)) {
    stop("Number of clusters 'ncl' must be a single integer.")
  }
  ncl <- as.integer(ncl)
  if (ncl < 1L || ncl > length(object$I_T_Y)) {
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