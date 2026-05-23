#' Compute py_x for new observations against training data
#'
#' Evaluates the kernel of new observations against a fixed training set
#' using the same bandwidth and kernel choices as the original fit. Used
#' by \code{predict.gibclust} to compute the conditional distribution
#' \eqn{p(y \mid x_{\text{new}})} required for assigning new observations
#' to clusters.
#'
#' @param X_train Training data frame (the data used to fit the model).
#' @param X_new Data frame of new observations to be evaluated.
#' @param s Bandwidth parameter(s) for continuous variables.
#' @param cat_cols Indices of categorical columns.
#' @param cont_cols Indices of continuous columns.
#' @param lambda Bandwidth parameter(s) for categorical variables.
#' @param contkernel Continuous kernel.
#' @param nomkernel Nominal kernel.
#' @param ordkernel Ordinal kernel.
#'
#' @return A list with components:
#'   \item{py_x_new}{An \eqn{n_{\text{train}} \times n_{\text{new}}} matrix
#'     of normalised kernel weights, where column \eqn{j} gives the
#'     distribution over training observations for new observation \eqn{j}.}
#'
#' @keywords internal
#' @noRd
coord_to_pxy_eval_R <- function(X_train, X_new, s, cat_cols, cont_cols, lambda,
                                contkernel = "gaussian",
                                nomkernel = "aitchisonaitken",
                                ordkernel = "liracine") {
  bws <- rep(NA, ncol(X_train))
  bws[cat_cols]  <- lambda
  bws[cont_cols] <- s
  py_x_new <- np::npksum(bws = bws,
                         txdat = X_train,
                         exdat = X_new,
                         ckertype = contkernel,
                         ukertype = nomkernel,
                         okertype = ordkernel,
                         return.kernel.weights = TRUE)$kw
  col_sums <- colSums(py_x_new)
  # Warning for new obs with zero kernel density to all training data
  zero_cols <- which(col_sums == 0)
  if (length(zero_cols) > 0) {
    warning(sprintf("Observation(s) %s in newdata have zero kernel density against the training set. Assignments for these may be unreliable.",
                    paste(zero_cols, collapse = ", ")))
    col_sums[zero_cols] <- 1
    py_x_new[, zero_cols] <- 1 / nrow(py_x_new)
  }
  py_x_new <- sweep(py_x_new, 2, col_sums, "/")
  
  list(py_x_new = py_x_new)
}

#' Reconstruct q(y | t) from training data and cluster assignments
#'
#' Given the training conditional distribution \eqn{p(y \mid x)} and a fit's
#' cluster assignments, computes the cluster profile distribution
#' \eqn{q(y \mid t)} for each cluster \eqn{t}. For hard fits, this is the
#' cluster-wise weighted average of training columns; for soft fits, it uses
#' the membership matrix as weights.
#'
#' @param py_x The training \eqn{n_{\text{train}} \times n_{\text{train}}}
#'   conditional distribution matrix.
#' @param cluster Either an integer vector of length \eqn{n_{\text{train}}}
#'   (hard) or an \eqn{n_{\text{cl}} \times n_{\text{train}}} membership
#'   matrix (soft).
#' @param px The observation marginal (length \eqn{n_{\text{train}}}).
#' @param ncl Number of clusters.
#'
#' @return An \eqn{n_{\text{train}} \times n_{\text{cl}}} matrix whose
#'   columns are the cluster profiles \eqn{q(y \mid t)}.
#'
#' @keywords internal
#' @noRd
.reconstruct_qy_t <- function(py_x, cluster, px, ncl) {
  n_train <- ncol(py_x)
  px_vec <- as.numeric(px)
  if (is.matrix(cluster)) {
    qt_x <- cluster
    if (nrow(qt_x) != ncl || ncol(qt_x) != n_train) {
      stop("Membership matrix dimensions do not match (ncl, n_train).")
    }
  } else {
    qt_x <- matrix(0, nrow = ncl, ncol = n_train)
    for (k in seq_len(ncl)) {
      qt_x[k, cluster == k] <- 1
    }
  }
  qt <- as.numeric(qt_x %*% px_vec)
  qy_t_unnorm <- py_x %*% (px_vec * t(qt_x))
  qy_t <- sweep(qy_t_unnorm, 2, qt, "/")
  empty <- which(qt == 0)
  if (length(empty) > 0) {
    warning(sprintf("Cluster(s) %s are empty; using uniform fallback.",
                    paste(empty, collapse = ", ")))
    qy_t[, empty] <- 1 / n_train
  }
  
  qy_t
}