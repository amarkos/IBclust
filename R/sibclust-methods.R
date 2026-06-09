#' Methods for sibclust objects
#'
#' S3 methods available for \code{"sibclust"} objects, including extractors
#' for the cluster assignments and model parameters, an information-metrics
#' accessor, a prediction method for new data, and diagnostic plotting.
#'
#' @details
#' The following methods are available:
#' \itemize{
#'   \item \code{\link[=print.sibclust]{print}} and \code{\link[=summary.sibclust]{summary}}:
#'     concise and detailed descriptions of the cluster solution.
#'   \item \code{\link[=fitted.sibclust]{fitted}}: extract cluster assignments.
#'     The argument \code{method = "classes"} (default) returns hard cluster
#'     labels; \code{method = "soft"} returns the equivalent one-hot
#'     membership matrix.
#'   \item \code{\link[=coef.sibclust]{coef}}: extract the model's bandwidth
#'     hyperparameters (\code{s}, \code{lambda}).
#'   \item \code{\link[=info_metrics]{info_metrics}}: extract
#'     information-theoretic quantities \eqn{H(T)}, \eqn{H(T \mid X)},
#'     \eqn{I(T; X)}, and \eqn{I(T; Y)}.
#'   \item \code{\link[=predict.sibclust]{predict}}: assign new observations to
#'     clusters via the sequential merge-cost (Jensen--Shannon) rule.
#'   \item \code{\link[=plot.sibclust]{plot}}: produce diagnostic plots
#'     (\code{type = "sizes"}, \code{"info"}, \code{"importance"}, or
#'     \code{"similarity"}).
#' }
#'
#' @name sibclust-methods
#' @aliases print.sibclust summary.sibclust print.summary.sibclust plot.sibclust fitted.sibclust coef.sibclust
#' @seealso \code{\link{sIBmix}}
#' @importFrom graphics barplot points
NULL

#' @rdname sibclust-methods
#' @param x a \code{sibclust} object.
#' @param object a \code{sibclust} object.
#' @param ... additional arguments passed to or from other methods.
#' @method print sibclust
#' @exportS3Method
print.sibclust <- function(x, ...) {
  header <- "Hard clustering with sIBmix"
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n")
  
  cat("Call: "); print(x$call)
  cat(sprintf("Observations: %d   Clusters: %d\n", x$n, x$ncl))
  cat(sprintf("Continuous variables: %d   Categorical variables: %d\n",
              length(x$contcols), length(x$catcols)))
  
  mi <- x$MutualInfo; ent <- x$Entropy; cent <- x$CondEntropy
  cat(sprintf("Mutual information I(Y;T): %s\n",
              if (is.finite(mi)) sprintf("%.6f", mi) else "Inf"))
  cat(sprintf("Entropy H(T): %s Conditional entropy H(T|X): %s\n",
              if (is.finite(ent)) sprintf("%.6f", ent) else "NA",
              if (is.finite(cent)) sprintf("%.6f", cent) else "NA"))
  
  cat(sprintf("Converged: %s\n", if (isTRUE(x$converged)) "TRUE" else "FALSE"))
  if (!is.na(x$iters)) cat(sprintf("Sweeps: %d\n", x$iters))
  if (!is.na(x$eps))   cat(sprintf("Convergence tolerance (eps): %g\n", x$eps))
  
  cl <- as.integer(x$Cluster)
  if (length(cl) == x$n) {
    cat("Cluster sizes:\n")
    print(sort(table(cl), decreasing = TRUE))
  } else {
    cat("Cluster labels unavailable or malformed.\n")
  }
  invisible(x)
}

#' @rdname sibclust-methods
#' @method summary sibclust
#' @exportS3Method
summary.sibclust <- function(object, ...) {
  cl <- as.integer(object$Cluster)
  sizes <- if (length(cl) == object$n) sort(table(cl), decreasing = TRUE) else integer(0)
  
  out <- list(
    call = object$call,
    n = object$n,
    ncl = object$ncl,
    sizes = sizes,
    proportions = if (length(sizes)) round(sizes / sum(sizes), 4) else numeric(0),
    Entropy = object$Entropy,
    CondEntropy = object$CondEntropy,
    MutualInfo = object$MutualInfo,
    s = object$s,
    lambda = object$lambda,
    contcols = object$contcols,
    catcols = object$catcols,
    kernels = object$kernels,
    iters = object$iters,
    converged = isTRUE(object$converged),
    eps = object$eps
  )
  class(out) <- "summary.sibclust"
  out
}

#' @rdname sibclust-methods
#' @method print summary.sibclust
#' @exportS3Method
print.summary.sibclust <- function(x, ...) {
  header <- "Summary of sIBmix clustering"
  cat(header, "\n", strrep("-", nchar(header)), "\n", sep = "")
  cat("Call: "); print(x$call)
  cat(sprintf("n = %d, k = %d\n\n", x$n, x$ncl))
  cat(sprintf("Continuous variables: %d   Categorical variables: %d\n\n",
              length(x$contcols), length(x$catcols)))
  
  if (length(x$sizes)) {
    cat("Cluster sizes:\n"); print(x$sizes)
    cat("\nProportions:\n"); print(x$proportions)
  } else {
    cat("Cluster labels unavailable.\n")
  }
  
  cat("\nInformation metrics:\n")
  cat(sprintf("Entropy H(T): %s\n",
              if (is.finite(x$Entropy)) sprintf("%.6f", x$Entropy) else "NA"))
  cat(sprintf("Conditional H(T|X): %s\n",
              if (is.finite(x$CondEntropy)) sprintf("%.6f", x$CondEntropy) else "NA"))
  cat(sprintf("Mutual Information I(Y;T): %s\n",
              if (is.finite(x$MutualInfo)) sprintf("%.6f", x$MutualInfo) else "Inf"))
  
  .format_vec <- function(v, name, digits = 4, max_show = 6) {
    if (!length(v)) return()
    r <- round(v, digits)
    if (length(r) > max_show) {
      cat(sprintf("%s = %s, ... (%d total)\n", name,
                  paste(r[1:max_show], collapse = ", "), length(r)))
    } else {
      cat(sprintf("%s = %s\n", name, paste(r, collapse = ", ")))
    }
  }
  
  cat("\nHyperparameters & details:\n")
  .format_vec(x$s, "s")
  .format_vec(x$lambda, "lambda")
  ks <- x$kernels
  cat(sprintf("Kernels = cont:%s, nom:%s, ord:%s\n", ks$cont, ks$nom, ks$ord))
  
  cat(sprintf("\nConverged: %s\n", if (isTRUE(x$converged)) "TRUE" else "FALSE"))
  if (!is.na(x$iters)) cat(sprintf("Sweeps: %d\n", x$iters))
  if (!is.na(x$eps))   cat(sprintf("Convergence tolerance (eps): %g\n", x$eps))
  invisible(x)
}

#' @rdname sibclust-methods
#' @param type Plot type: \code{"sizes"} (cluster sizes), \code{"info"} (information metrics), \code{"importance"} (variable importance bar chart), or \code{"similarity"} (heatmap of the kernel similarity matrix \eqn{P_{Y|X}}).
#' @param main Optional title.
#' @param col Optional color (or, for \code{type = "similarity"}, a colour palette vector).
#' @param X Original data frame used to fit \code{x}; required for
#'   \code{type = "importance"} and \code{type = "similarity"} when the
#'   fit was constructed with \code{keep_data = FALSE}. If the fit already
#'   contains the training data (\code{keep_data = TRUE}), any supplied
#'   \code{X} is ignored with a warning.
#' @param color_by_type Logical; if \code{TRUE}, colour bars by variable type
#'   (continuous / nominal / ordinal). Defaults to \code{TRUE}. Used only by
#'   \code{type = "importance"}.
#' @param order_by_cluster Logical; if \code{TRUE} (default), rows and columns
#'   of the similarity matrix are reordered by cluster assignment and
#'   cluster-boundary boxes are drawn. Used only by \code{type = "similarity"}.
#' @param colorbar Logical; if \code{TRUE} (default), draw a horizontal colour
#'   scale below the similarity heatmap. Used only by \code{type = "similarity"}.
#' @method plot sibclust
#' @exportS3Method
plot.sibclust <- function(x, type = c("sizes", "info", "importance", "similarity"),
                          X = NULL, color_by_type = TRUE, col = NULL,
                          order_by_cluster = TRUE, colorbar = TRUE,
                          main = NULL, ...) {
  type <- match.arg(type)
  
  if (type == "importance") {
    if (!is.null(x$training_data)) {
      if (!is.null(X)) {
        warning("Argument 'X' was supplied but the fitted object already contains the training data (keep_data = TRUE). Using the stored training data; the supplied 'X' will be ignored.")
      }
      X <- x$training_data
    } else if (is.null(X)) {
      stop("Argument 'X' (the original data frame) is required for type = 'importance', or refit with keep_data = TRUE.")
    }
    if (nrow(X) != x$n) {
      stop(sprintf("nrow(X) = %d does not match the fitted model's n = %d.",
                   nrow(X), x$n))
    }
    iyt <- .compute_variable_importance(
      X = X,
      cluster = as.integer(x$Cluster),
      s = x$s,
      lambda = x$lambda,
      contcols = x$contcols,
      catcols = x$catcols,
      kernels = x$kernels,
      nystrom_landmarks = NULL,
      scale = x$scale
    )
    .plot_variable_importance(iyt, X = X,
                              color_by_type = color_by_type,
                              col = col,
                              main = main, ...)
    return(invisible(x))
  }
  
  if (type == "similarity") {
    if (!is.null(x$training_data)) {
      if (!is.null(X)) {
        warning("Argument 'X' was supplied but the fitted object already contains the training data (keep_data = TRUE). Using the stored training data; the supplied 'X' will be ignored.")
      }
      X <- x$training_data
    } else if (is.null(X)) {
      stop("Argument 'X' (the original data frame) is required for type = 'similarity', or refit with keep_data = TRUE.")
    }
    if (nrow(X) != x$n) {
      stop(sprintf("nrow(X) = %d does not match the fitted model's n = %d.",
                   nrow(X), x$n))
    }
    contcols <- x$contcols
    catcols <- x$catcols
    X <- as.data.frame(X)
    if (length(contcols) > 0 && isTRUE(x$scale)) {
      X[, contcols] <- scale(X[, contcols, drop = FALSE])
    }
    if (length(catcols) > 0) {
      X[, catcols] <- preprocess_cat_data(X[, catcols, drop = FALSE])
    }
    pxy_list <- coord_to_pxy_R(
      X = X,
      s = if (length(contcols) > 0L) x$s else -1,
      lambda = if (length(catcols)  > 0L) x$lambda else -1,
      cat_cols = catcols,
      cont_cols = contcols,
      contkernel = x$kernels$cont,
      nomkernel = x$kernels$nom,
      ordkernel = x$kernels$ord
    )
    M <- pxy_list$py_x
    cluster_sizes <- NULL
    if (isTRUE(order_by_cluster)) {
      cl <- as.integer(x$Cluster)
      ord <- order(cl)
      M <- M[ord, ord, drop = FALSE]
      cluster_sizes <- as.integer(table(cl))
    }
    zmax <- quantile(M[M > 0], 0.99, na.rm = TRUE)
    plot_col <- if (is.null(col)) colorRampPalette(c("white", "steelblue"))(100) else col
    if (is.null(main)) main <- expression("Similarity matrix " * P[Y * "|" * X])
    
    if (isTRUE(colorbar)) {
      old_mar <- par("mar")
      on.exit(par(mar = old_mar), add = TRUE)
      par(mar = c(5, 2, 4, 2) + 0.1)
    }
    image(t(M)[, nrow(M):1],
          col = plot_col,
          zlim = c(0, zmax),
          main = main,
          xlab = "", ylab = "", axes = FALSE,
          ...)
    box(lwd = 1, col = "black")
    if (!is.null(cluster_sizes)) {
      n <- nrow(M)
      ends <- cumsum(cluster_sizes)
      starts <- c(1, head(ends, -1) + 1)
      for (k in seq_along(cluster_sizes)) {
        x0 <- (starts[k] - 1) / n
        x1 <- ends[k] / n
        y0 <- 1 - ends[k] / n
        y1 <- 1 - (starts[k] - 1) / n
        rect(x0, y0, x1, y1, border = "black", lwd = 2.5)
      }
    }
    if (isTRUE(colorbar)) {
      usr <- par("usr")
      bar_x <- seq(usr[1], usr[2], length.out = length(plot_col) + 1)
      bar_y_top <- usr[3] - 0.04 * (usr[4] - usr[3])
      bar_y_bot <- usr[3] - 0.10 * (usr[4] - usr[3])
      for (i in seq_along(plot_col)) {
        rect(bar_x[i], bar_y_bot, bar_x[i + 1], bar_y_top,
             col = plot_col[i], border = NA, xpd = NA)
      }
      rect(usr[1], bar_y_bot, usr[2], bar_y_top,
           border = "black", lwd = 1, xpd = NA)
      tick_vals <- pretty(c(0, zmax), n = 5)
      tick_vals <- tick_vals[tick_vals > 0 & tick_vals < zmax]
      tick_vals <- c(0, tick_vals, zmax)
      tick_x <- usr[1] + (tick_vals / zmax) * (usr[2] - usr[1])
      text(tick_x, bar_y_bot - 0.025 * (usr[4] - usr[3]),
           labels = formatC(tick_vals, format = "g", digits = 2),
           cex = 0.7, xpd = NA, adj = c(0.5, 1))
      segments(tick_x, bar_y_bot, tick_x, bar_y_bot - 0.01 * (usr[4] - usr[3]),
               xpd = NA)
    }
    return(invisible(x))
  }
  
  if (type == "sizes") {
    cl <- as.integer(x$Cluster)
    if (length(cl) == x$n) {
      if (is.null(main)) main <- "Cluster sizes (sIBmix)"
      barplot(table(cl), ylab = "Count", xlab = "Cluster", main = main,
              col = if (is.null(col)) "gray" else col, ...)
    } else {
      warning("Hard labels unavailable or malformed.")
    }
  } else { # type == "info"
    vals <- c(`H(T)`   = x$Entropy,
              `H(T|X)` = x$CondEntropy,
              `I(Y;T)` = x$MutualInfo)
    vals[!is.finite(vals)] <- NA_real_
    if (is.null(main)) main <- "Information summary"
    barplot(vals, ylab = "Value", main = main,
            col = if (is.null(col)) "gray" else col, ...)
  }
  invisible(x)
}

#' @rdname sibclust-methods
#' @param method For \code{fitted}: either \code{"classes"} (default; the
#'   integer cluster-label vector) or \code{"soft"} (the equivalent one-hot
#'   binary membership matrix).
#' @method fitted sibclust
#' @exportS3Method
fitted.sibclust <- function(object, method = c("classes", "soft"), ...) {
  method <- match.arg(method)
  cl <- as.integer(object$Cluster)
  if (method == "classes") return(cl)
  
  ncl <- object$ncl
  M <- matrix(0, nrow = ncl, ncol = length(cl))
  for (k in seq_len(ncl)) M[k, cl == k] <- 1
  M
}

#' @rdname sibclust-methods
#' @method coef sibclust
#' @exportS3Method
coef.sibclust <- function(object, ...) {
  out <- list()
  if (length(object$contcols) > 0L) out$s <- object$s
  if (length(object$catcols)  > 0L) out$lambda <- object$lambda
  out
}

#' Predict cluster assignments for new observations
#'
#' Assigns new observations to clusters using a fitted \code{sibclust} model.
#' Each new observation is assigned to the cluster minimising the sequential
#' merge cost, i.e. the (prior-weighted) loss of information incurred by
#' absorbing the observation into the cluster,
#' \eqn{(p_x + p_t)\,\mathrm{JS}_\pi(p(Y \mid x_{\text{new}}), q(Y \mid t))}, i.e.
#' the same criterion optimised during fitting.
#'
#' Continuous variables in \code{newdata} are scaled using the means and
#' standard deviations of the training data \code{X} when \code{object$scale}
#' is \code{TRUE}. Categorical variables are relabelled to match the encoding
#' used at fit time; new observations with categorical levels unseen in the
#' training data produce an error.
#'
#' @param object A fitted \code{sibclust} object produced by \code{\link{sIBmix}}.
#' @param newdata A data frame of new observations to be assigned. Must have
#'   the same columns (matching names and order) as the training data. If
#'   \code{NULL} (default), predictions are returned for the training data
#'   itself.
#' @param X The original training data frame. Optional if \code{object} was
#'   constructed with \code{keep_data = TRUE}; required otherwise.
#' @param ... Additional arguments (currently ignored).
#'
#' @return An integer vector of length \code{nrow(newdata)} giving the
#'   predicted cluster label for each new observation.
#'
#' @seealso \code{\link{sIBmix}}, \code{\link[stats]{fitted}}.
#'
#' @method predict sibclust
#' @exportS3Method
predict.sibclust <- function(object, newdata = NULL, X = NULL, ...) {
  if (is.null(X)) {
    if (is.null(object$training_data)) {
      stop("'X' (the training data) must be supplied. Alternatively, refit with keep_data = TRUE to store training data in the model object.")
    }
    X <- object$training_data
  }
  if (is.null(newdata)) newdata <- X
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  if (!is.data.frame(X)) X <- as.data.frame(X)
  
  if (!identical(names(newdata), names(X))) {
    stop("Columns of 'newdata' must match columns of 'X' exactly (same names, same order).")
  }
  if (nrow(X) != object$n) {
    stop(sprintf("nrow(X) = %d does not match the fitted model's n = %d.",
                 nrow(X), object$n))
  }
  
  contcols <- object$contcols
  catcols  <- object$catcols
  
  if (length(contcols) > 0L && isTRUE(object$scale)) {
    train_means <- colMeans(X[, contcols, drop = FALSE])
    train_sds   <- apply(X[, contcols, drop = FALSE], 2, stats::sd)
    X[, contcols] <- scale(X[, contcols, drop = FALSE],
                           center = train_means, scale = train_sds)
    newdata[, contcols] <- scale(newdata[, contcols, drop = FALSE],
                                 center = train_means, scale = train_sds)
  }
  if (length(catcols) > 0L) {
    for (j in catcols) {
      lev <- levels(factor(X[[j]]))
      X[[j]]       <- as.integer(factor(X[[j]], levels = lev))
      newdata[[j]] <- as.integer(factor(newdata[[j]], levels = lev))
    }
  }
  
  eval_list <- coord_to_pxy_eval_R(
    X_train = X, X_new = newdata,
    s = if (length(contcols) > 0L) object$s else -1,
    lambda = if (length(catcols)  > 0L) object$lambda else -1,
    cat_cols = catcols, cont_cols = contcols,
    contkernel = object$kernels$cont,
    nomkernel = object$kernels$nom,
    ordkernel = object$kernels$ord
  )
  py_x_new <- eval_list$py_x_new
  
  pxy_train <- coord_to_pxy_R(
    X = X,
    s = if (length(contcols) > 0L) object$s else -1,
    lambda = if (length(catcols)  > 0L) object$lambda else -1,
    cat_cols = catcols, cont_cols = contcols,
    contkernel = object$kernels$cont,
    nomkernel = object$kernels$nom,
    ordkernel = object$kernels$ord
  )
  py_x_train <- pxy_train$py_x
  px_train <- pxy_train$px
  
  ncl <- object$ncl
  qy_t <- .reconstruct_qy_t(py_x_train, object$Cluster, px_train, ncl)
  
  # cluster masses and single-point mass
  qt <- as.numeric(table(factor(object$Cluster, levels = seq_len(ncl))) / object$n)
  p_x <- 1 / object$n
  
  H <- function(p) entropySingle(pmax(as.numeric(p), 0))
  Hc <- vapply(seq_len(ncl), function(t) H(qy_t[, t]), numeric(1))
  
  n_new <- ncol(py_x_new)
  score <- matrix(0, nrow = ncl, ncol = n_new)
  for (t in seq_len(ncl)) {
    pt <- qt[t]
    denom <- p_x + pt
    for (j in seq_len(n_new)) {
      m_t <- (p_x * py_x_new[, j] + pt * qy_t[, t]) / denom
      # exact merge cost up to j-constant term p_x H(p(Y|x_new)):
      score[t, j] <- denom * H(m_t) - pt * Hc[t]
    }
  }
  apply(score, 2, which.min)
}