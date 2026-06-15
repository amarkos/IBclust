#' Methods for gibclust objects
#'
#' S3 methods available for \code{gibclust} objects, including extractors
#' for the cluster assignments and model parameters, an information-metrics
#' accessor, a prediction method for new data, and diagnostic plotting.
#'
#' @details
#' The following methods are available:
#' \itemize{
#'   \item \code{\link[=print.gibclust]{print}} and \code{\link[=summary.gibclust]{summary}}:
#'     concise and detailed descriptions of the cluster solution.
#'   \item \code{\link[=fitted.gibclust]{fitted}}: extract cluster assignments.
#'     The argument \code{method = "classes"} (default) returns hard cluster
#'     labels; \code{method = "soft"} returns the membership matrix (one-hot
#'     for DIBmix; fuzzy for IBmix and GIBmix).
#'   \item \code{\link[=coef.gibclust]{coef}}: extract the model's
#'     hyperparameters as a list, including bandwidths (\code{s},
#'     \code{lambda}) and IB parameters (\code{beta}, \code{alpha}).
#'   \item \code{\link[=info_metrics]{info_metrics}}: extract
#'     information-theoretic quantities \eqn{H(T)}, \eqn{H(T \mid X)},
#'     \eqn{I(T; X)}, and \eqn{I(T; Y)}.
#'   \item \code{\link[=predict.gibclust]{predict}}: assign new observations to
#'     clusters. For DIBmix fits, returns hard labels via the minimisation of
#'     Kullback--Leibler divergence; for IBmix and GIBmix fits, returns a
#'     soft membership matrix via the soft IB assignment rule.
#'   \item \code{\link[=plot.gibclust]{plot}}: produce diagnostic plots
#'     (\code{type = "sizes"}, \code{"info"}, \code{"beta"},
#'     \code{"importance"}, \code{"similarity"}, or \code{"membership"}).
#' }
#'
#' @name gibclust-methods
#' @aliases print.gibclust summary.gibclust print.summary.gibclust plot.gibclust fitted.gibclust coef.gibclust
#' @seealso \code{\link{DIBmix}}, \code{\link{IBmix}}, \code{\link{GIBmix}}
#' @importFrom graphics barplot points
NULL

#' @rdname gibclust-methods
#' @param x a \code{gibclust} object.
#' @param object a \code{gibclust} object.
#' @param ... additional arguments passed to or from other methods.
#' @method print gibclust
#' @exportS3Method
print.gibclust <- function(x, ...) {
  # Header depends on alpha
  header <- if (isTRUE(all.equal(x$alpha, 1))) {
    "Fuzzy clustering with IBmix"
  } else if (isTRUE(all.equal(x$alpha, 0))) {
    "Hard clustering with DIBmix"
  } else {
    "Fuzzy clustering with GIBmix"
  }
  cat(header, "\n")
  cat(strrep("-", nchar(header)), "\n")
  
  cat("Call: "); print(x$call)
  cat(sprintf("Observations: %d   Clusters: %d\n", x$n, x$ncl))
  cat(sprintf("Continuous variables: %d   Categorical variables: %d\n",
              length(x$contcols), length(x$catcols)))
  
  mi <- x$MutualInfo
  ent <- x$Entropy
  cent <- x$CondEntropy
  
  cat(sprintf("Mutual information I(Y;T): %s\n",
              if (is.finite(mi)) sprintf("%.6f", mi) else "Inf"))
  cat(sprintf("Entropy H(T): %s Conditional entropy H(T|X): %s\n",
              if (is.finite(ent)) sprintf("%.6f", ent) else "NA",
              if (is.finite(cent)) sprintf("%.6f", cent) else "NA"))
  
  cat(sprintf("Converged: %s\n", if (isTRUE(x$converged)) "TRUE" else "FALSE"))
  if (!is.na(x$iters)) cat(sprintf("Iterations: %d\n", x$iters))
  if (!is.na(x$conv_tol)) cat(sprintf("Convergence tolerance: %g\n", x$conv_tol))
  
  # Harden fuzzy memberships: C is ncl x n (rows = clusters, cols = obs)
  .harden_membership <- function(C) apply(C, 2, which.max)
  
  if (isTRUE(all.equal(x$alpha, 0))) {
    # DIBmix: hard labels vector expected
    cl <- x$Cluster
    if (is.integer(cl) || is.numeric(cl)) cl <- as.integer(cl)
    if (is.vector(cl) && length(cl) == x$n) {
      tab <- sort(table(cl), decreasing = TRUE)
      cat("Cluster sizes:\n")
      print(tab)
    } else {
      cat("Cluster labels unavailable or malformed.\n")
    }
  } else {
    # IB/GIB: fuzzy membership matrix expected
    C <- x$Cluster
    if (is.matrix(C) && nrow(C) == x$ncl && ncol(C) == x$n) {
      # Mean membership per cluster
      avg_mem <- rowMeans(C)
      hard_lab <- .harden_membership(C)
      hard_sizes <- sort(table(hard_lab), decreasing = TRUE)
      
      cat("Fuzzy memberships:\n")
      cat("Mean membership per cluster:\n")
      print(round(avg_mem, 4))
      cat("Hardened sizes (argmax):\n")
      print(hard_sizes)
    } else {
      cat("Membership matrix unavailable or malformed.\n")
    }
  }
  invisible(x)
}

#' @rdname gibclust-methods
#' @method summary gibclust
#' @exportS3Method
summary.gibclust <- function(object, ...) {
  variant <- if (isTRUE(all.equal(object$alpha, 0))) {
    "DIBmix"
  } else if (isTRUE(all.equal(object$alpha, 1))) {
    "IBmix"
  } else {
    "GIBmix"
  }
  
  # helpers for fuzzy case (C is ncl x n)
  harden <- function(C) apply(C, 2, which.max)
  
  if (variant == "DIBmix") {
    cl <- object$Cluster
    sizes <- if (is.vector(cl) && length(cl) == object$n) {
      sort(table(cl), decreasing = TRUE)
    } else integer(0)
    
    out <- list(
      call = object$call,
      n = object$n,
      ncl = object$ncl,
      variant = variant,
      sizes = sizes,
      proportions = if (length(sizes)) round(sizes / sum(sizes), 4) else numeric(0),
      Entropy = object$Entropy,
      CondEntropy = object$CondEntropy,
      MutualInfo = object$MutualInfo,
      beta = object$beta,
      alpha = object$alpha,
      s = object$s,
      lambda = object$lambda,
      contcols = object$contcols,
      catcols = object$catcols,
      kernels = object$kernels,
      iters = object$iters,
      converged = isTRUE(object$converged),
      conv_tol = object$conv_tol
    )
  } else {
    C <- object$Cluster
    mean_membership <- if (is.matrix(C) && nrow(C) == object$ncl && ncol(C) == object$n) {
      rowMeans(C)
    } else rep(NA_real_, object$ncl)
    
    hardened <- if (is.matrix(C)) harden(C) else rep(NA_integer_, object$n)
    hardened_sizes <- if (all(is.finite(hardened))) sort(table(hardened), decreasing = TRUE) else integer(0)
    
    out <- list(
      call = object$call,
      n = object$n,
      ncl = object$ncl,
      variant = variant,
      mean_membership = mean_membership,
      hardened_sizes = hardened_sizes,
      Entropy = object$Entropy,
      CondEntropy = object$CondEntropy,
      MutualInfo = object$MutualInfo,
      beta = object$beta,
      alpha = object$alpha,
      s = object$s,
      lambda = object$lambda,
      contcols = object$contcols,
      catcols = object$catcols,
      kernels = object$kernels,
      iters = object$iters,
      converged = isTRUE(object$converged),
      conv_tol = object$conv_tol
    )
  }
  
  class(out) <- "summary.gibclust"
  out
}

#' @rdname gibclust-methods
#' @method print summary.gibclust
#' @exportS3Method
print.summary.gibclust <- function(x, ...) {
  header <- switch(x$variant,
                   "DIBmix" = "Summary of DIBmix clustering",
                   "IBmix"  = "Summary of IBmix clustering",
                   "GIBmix" = "Summary of GIBmix clustering",
                   "Summary of GIBmix clustering")
  cat(header, "\n", strrep("-", nchar(header)), "\n", sep = "")
  cat("Call: ")
  print(x$call)
  cat(sprintf("n = %d, k = %d\n\n", x$n, x$ncl))
  cat(sprintf("Continuous variables: %d   Categorical variables: %d\n\n",
              length(x$contcols), length(x$catcols)))
  
  # Cluster summaries
  if (identical(x$variant, "DIBmix")) {
    if (length(x$sizes)) {
      cat("Cluster sizes:\n"); print(x$sizes)
      cat("\nProportions:\n"); print(x$proportions)
    } else {
      cat("Cluster labels unavailable.\n")
    }
  } else {
    if (length(x$mean_membership)) {
      cat("Mean membership per cluster:\n")
      print(round(x$mean_membership, 4))
    }
    if (length(x$hardened_sizes)) {
      cat("\nHardened sizes (argmax):\n")
      print(x$hardened_sizes)
    }
  }
  
  # Info metrics
  cat("\nInformation metrics:\n")
  cat(sprintf("Entropy H(T): %s\n",
              if (is.finite(x$Entropy)) sprintf("%.6f", x$Entropy) else "NA"))
  cat(sprintf("Conditional H(T|X): %s\n",
              if (is.finite(x$CondEntropy)) sprintf("%.6f", x$CondEntropy) else "NA"))
  cat(sprintf("Mutual Information I(Y;T): %s\n",
              if (is.finite(x$MutualInfo)) sprintf("%.6f", x$MutualInfo) else "Inf"))
  
  # Compact vector printing helper
  .format_vec <- function(v, name, digits = 4, max_show = 6) {
    if (!length(v)) return()
    rounded <- round(v, digits)
    if (length(rounded) > max_show) {
      shown <- paste(rounded[1:max_show], collapse = ", ")
      cat(sprintf("%s = %s, ... (%d total)\n", name, shown, length(rounded)))
    } else {
      cat(sprintf("%s = %s\n", name, paste(rounded, collapse = ", ")))
    }
  }
  
  cat("\nHyperparameters & details:\n")
  .format_vec(x$beta, "beta")
  .format_vec(x$s, "s")
  .format_vec(x$lambda, "lambda")
  cat(sprintf("alpha = %s\n", deparse(x$alpha)))
  ks <- x$kernels
  cat(sprintf("Kernels = cont:%s, nom:%s, ord:%s\n", ks$cont, ks$nom, ks$ord))
  
  cat(sprintf("\nConverged: %s\n", if (isTRUE(x$converged)) "TRUE" else "FALSE"))
  if (!is.na(x$iters)) cat(sprintf("Iterations: %d\n", x$iters))
  if (!is.na(x$conv_tol)) cat(sprintf("Convergence tolerance: %g\n", x$conv_tol))
  invisible(x)
}

#' @rdname gibclust-methods
#' @param type Plot type: \code{"sizes"} (cluster sizes), \code{"info"} (information metrics), \code{"beta"} (log(beta) trajectory; DIBmix only), \code{"importance"} (variable importance bar chart), \code{"similarity"} (heatmap of the kernel similarity matrix \eqn{P_{Y|X}}), or \code{"membership"} (parallel coordinates plot illustrating cluster membership).
#' @param main Optional title.
#' @param col Optional color (or, for type = "similarity" and type = "membership", a colour palette vector).
#' @param X Original data frame used to fit \code{x}; required for
#'   \code{type = "importance"} and \code{type = "similarity"} when the
#'   fit was constructed with \code{keep_data = FALSE}. If the fit already
#'   contains the training data (\code{keep_data = TRUE}), any supplied
#'   \code{X} is ignored with a warning.
#' @param color_by_type Logical; if \code{TRUE}, colour bars by variable type
#'   (continuous / nominal / ordinal). Defaults to \code{TRUE}. Used only by
#'   \code{type = "importance"}.
#' @param order_by_cluster Logical; if \code{TRUE} (default), rows and columns
#'   of the similarity matrix are reordered by cluster assignment and cluster-boundary
#'   boxes are drawn. Used only by \code{type = "similarity"}.
#' @param groups Optional grouping vector of length \eqn{n} used to colour
#'   the lines in \code{type = "membership"} (e.g. an external class label).
#'   If \code{NULL} (default), lines are coloured by the dominant (argmax)
#'   cluster. Used only by \code{type = "membership"}.
#' @param colorbar Logical; if \code{TRUE} (default), draw a horizontal
#'   colour scale below the similarity heatmap. Used only by
#'   \code{type = "similarity"}.
#' @method plot gibclust
#' @exportS3Method
plot.gibclust <- function(x, type = c("sizes", "info", "beta", "importance", "similarity", "membership"),
                          X = NULL, color_by_type = TRUE, col = NULL,
                          order_by_cluster = TRUE, groups = NULL,
                          colorbar = TRUE,
                          main = NULL, ...) {
  type <- match.arg(type)
  harden <- function(C) apply(C, 2, which.max)
  
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
    cluster <- if (is.matrix(x$Cluster)) harden(x$Cluster) else x$Cluster
    iyt <- .compute_variable_importance(
      X = X,
      cluster = cluster,
      s = x$s,
      lambda = x$lambda,
      contcols = x$contcols,
      catcols = x$catcols,
      kernels = x$kernels,
      nystrom_landmarks = x$nystrom_landmarks,
      scale = x$scale
    )
    .plot_variable_importance(iyt, X = X,
                              color_by_type = color_by_type,
                              col = col,
                              main = main, ...)
    return(invisible(x))
  }
  
  if (type == "membership") {
    C <- x$Cluster
    if (is.matrix(C)) {
      U <- t(C)
    } else {
      warning("Membership plot is most informative for fuzzy fits (IBmix / GIBmix). For a hard (DIBmix) fit the memberships are binary.")
      cl_int <- as.integer(C)
      k <- if (!is.null(x$ncl)) x$ncl else max(cl_int)
      U <- matrix(0, nrow = length(cl_int), ncol = k)
      U[cbind(seq_along(cl_int), cl_int)] <- 1
    }
    k <- ncol(U)
    n <- nrow(U)
    dom <- max.col(U, ties.method = "first")
    
    if (!is.null(groups)) {
      if (length(groups) != n) {
        stop(sprintf("'groups' must have length n = %d.", n))
      }
      grp <- as.factor(groups)
    } else {
      grp <- factor(dom, levels = seq_len(k))
    }
    n_grp <- nlevels(grp)
    if (is.null(col)) {
      base_palette <- grDevices::hcl.colors(max(n_grp, 2), palette = "Dynamic")
      line_cols <- base_palette[as.integer(grp)]
    } else {
      line_cols <- rep_len(col, n_grp)[as.integer(grp)]
    }
    if (is.null(main)){
      main <- "Cluster membership probabilities"
    }
    plot(NA, xlim = c(1, k), ylim = c(0, 1),
         xlab = "Cluster", ylab = "Membership probability",
         main = main, xaxt = "n", ...)
    axis(1, at = seq_len(k), labels = seq_len(k))
    alpha_lines <- grDevices::adjustcolor(line_cols, alpha.f = 0.4)
    xs <- seq_len(k)
    for (i in seq_len(nrow(U))) {
      lines(xs, U[i, ], col = alpha_lines[i], lwd = 1)
    }
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
    if (!is.null(x$nystrom_landmarks)) {
      pxy_list <- coord_to_pxy_nystrom_R(
        X = X,
        s = if (length(contcols) > 0) x$s else -1,
        lambda = if (length(catcols)  > 0) x$lambda else -1,
        cat_cols = catcols,
        cont_cols = contcols,
        contkernel = x$kernels$cont,
        nomkernel = x$kernels$nom,
        ordkernel = x$kernels$ord,
        n_landmarks = length(x$nystrom_landmarks),
        landmark_indices = x$nystrom_landmarks
      )
      B <- pxy_list$py_x$B
      cs <- pxy_list$py_x$col_sums
      M <- tcrossprod(B)
      M <- sweep(M, 2, cs, FUN = "/")
    } else {
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
    }
    cluster_sizes <- NULL
    if (isTRUE(order_by_cluster)) {
      cl <- if (is.matrix(x$Cluster)) harden(x$Cluster) else as.integer(x$Cluster)
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
    if (isTRUE(all.equal(x$alpha, 0))) {
      cl <- x$Cluster
      if (is.vector(cl) && length(cl) == x$n) {
        tab <- table(as.integer(cl))
        if (is.null(main)) main <- "Cluster sizes (DIBmix)"
        barplot(tab, ylab = "Count", xlab = "Cluster", main = main,
                col = if (is.null(col)) "gray" else col, ...)
      } else {
        warning("Hard labels unavailable or malformed for DIBmix.")
      }
    } else {
      C <- x$Cluster
      if (is.matrix(C) && nrow(C) == x$ncl && ncol(C) == x$n) {
        lab <- harden(C)
        tab <- table(lab)
        if (is.null(main)) {
          main <- if (isTRUE(all.equal(x$alpha, 1))) "Hardened sizes (IBmix)"
          else "Hardened sizes (GIBmix)"
        }
        barplot(tab, ylab = "Count", xlab = "Cluster", main = main,
                col = if (is.null(col)) "gray" else col, ...)
      } else {
        warning("Membership matrix unavailable or malformed for IB/GIB.")
      }
    }
  } else if (type == "info") {
    vals <- c(`H(T)`   = x$Entropy,
              `H(T|X)` = x$CondEntropy,
              `I(Y;T)` = x$MutualInfo)
    vals[!is.finite(vals)] <- NA_real_
    if (is.null(main)) main <- "Information summary"
    barplot(vals, ylab = "Value", main = main,
            col = if (is.null(col)) "gray" else col, ...)
    
  } else {
    if (!isTRUE(all.equal(x$alpha, 0))) {
      warning("type='beta' is available only for DIBmix (alpha = 0).")
      return(invisible(x))
    }
    b <- x$beta
    if (!length(b)) {
      warning("No beta trajectory available to plot.")
      return(invisible(x))
    }
    idx <- 0:(length(b) - 1)
    logb <- log(b)
    line_col <- if (is.null(col)) "black" else col
    if (is.null(main)) main <- expression(log(beta) ~ " trajectory (DIBmix)")
    if (!all(is.finite(logb))) {
      warning("Non-finite values in log(beta); some points omitted.")
    }
    plot(idx, logb, type = "l", col = line_col,
         xlab = "Iteration", ylab = expression(log(beta)),
         main = main, ...)
    points(idx, logb, col = line_col, ...)
  }
  invisible(x)
}

#' @rdname gibclust-methods
#' @param method For \code{fitted}: either \code{"classes"}
#'   (default; hard cluster labels obtained via argmax of the membership
#'   matrix) or \code{"soft"} (the raw fuzzy membership matrix). For hard
#'   fits (DIBmix), \code{"classes"} returns the integer label vector and
#'   \code{"soft"} returns the equivalent one-hot binary matrix.
#' @method fitted gibclust
#' @exportS3Method
fitted.gibclust <- function(object, method = c("classes", "soft"), ...) {
  method <- match.arg(method)
  cl <- object$Cluster
  
  if (is.matrix(cl)) {
    if (method == "soft") return(cl)
    return(apply(cl, 2, which.max))
  }
  
  if (method == "classes") return(as.integer(cl))
  
  cl <- as.integer(cl)
  ncl <- object$ncl
  M <- matrix(0, nrow = ncl, ncol = length(cl))
  for (k in seq_len(ncl)) {
    M[k, cl == k] <- 1
  }
  M
}

#' @rdname gibclust-methods
#' @method coef gibclust
#' @exportS3Method
coef.gibclust <- function(object, ...) {
  out <- list()
  if (length(object$contcols) > 0L) out$s <- object$s
  if (length(object$catcols)  > 0L) out$lambda <- object$lambda
  out$beta  <- object$beta
  out$alpha <- object$alpha
  out
}

#' Predict cluster assignments for new observations
#'
#' Assigns new observations to clusters using a fitted \code{gibclust} model.
#' For hard fits (\code{DIBmix}, \code{alpha = 0}), returns integer cluster
#' labels based on the minimisation of the Kullback--Leibler divergence between the new
#' observation's conditional distribution and each cluster conditional density. For
#' soft fits (\code{IBmix}, \code{GIBmix}), returns the full membership matrix
#' based on the (G)IB update, with the fitted \eqn{\beta}.
#'
#' Prediction consists of three steps. First, the conditional distribution
#' \eqn{p(y \mid x_{\text{new}})} of each new observation over the training
#' data is computed using the same kernel and bandwidths as the fit. Second,
#' the cluster profiles \eqn{q(y \mid t)} are reconstructed from the training
#' data and the stored cluster assignments. Third, the Kullback--Leibler
#' divergence \eqn{D_{KL}(p(y \mid x_{\text{new}}) \,||\, q(y \mid t))} is
#' computed for each (new observation, cluster) pair. New observations are
#' then assigned via argmin (DIBmix) or via the soft assignment rule
#' \eqn{q(t \mid x_{\text{new}}) \propto q(t) \exp(-\beta \, D_{KL})}
#' (IBmix and GIBmix). 
#'
#' Continuous variables in \code{newdata} are scaled using the means and
#' standard deviations of the training data \code{X} when \code{object$scale}
#' is \code{TRUE}. Categorical variables are relabelled to match the encoding
#' used at fit time. New observations containing categorical levels not seen
#' in the training data will produce an error.
#'
#' @param object A fitted \code{gibclust} object produced by \code{\link{DIBmix}},
#'   \code{\link{IBmix}}, or \code{\link{GIBmix}}.
#' @param newdata A data frame of new observations to be assigned. Must have
#'   the same columns (with matching names and order) as the training data.
#'   If \code{NULL} (default), predictions are returned for the training data
#'   itself, providing a self-consistency check against \code{fitted(object)}.
#' @param X The original training data frame used to fit \code{object}.
#'   Optional if \code{object} was constructed with \code{keep_data = TRUE};
#'   in that case the stored training data is used automatically. Required
#'   otherwise.
#' @param ... Additional arguments (currently ignored).
#'
#' @return For DIBmix fits, an integer vector of length \code{nrow(newdata)}
#'   giving the predicted cluster label for each new observation. For IBmix
#'   and GIBmix fits, a numeric matrix of dimension \code{object$ncl} by
#'   \code{nrow(newdata)} containing soft membership probabilities, with
#'   columns summing to 1.
#'
#' @seealso \code{\link{DIBmix}}, \code{\link{IBmix}}, \code{\link{GIBmix}},
#'   \code{\link[stats]{fitted}}.
#'
#' @examples
#' # Hold-out prediction on the iris data
#' set.seed(1)
#' test_idx  <- sample(seq_len(nrow(iris)), 30)
#' train_idx <- setdiff(seq_len(nrow(iris)), test_idx)
#'
#' fit <- DIBmix(iris[train_idx, -5], ncl = 3, nstart = 5, keep_data = TRUE)
#' preds <- predict(fit, newdata = iris[test_idx, -5])
#' table(preds, iris$Species[test_idx])
#'
#' # Grid prediction for visualising the decision boundary
#' fit2 <- DIBmix(faithful, ncl = 3, nstart = 5, keep_data = TRUE)
#' grid <- expand.grid(
#'   eruptions = seq(min(faithful$eruptions), max(faithful$eruptions), length = 50),
#'   waiting = seq(min(faithful$waiting), max(faithful$waiting), length = 50)
#' )
#' grid_pred <- predict(fit2, newdata = grid)
#' plot(grid, col = grid_pred + 1, pch = 15, cex = 0.6)
#' points(faithful, pch = fit2$Cluster, col = "black")
#'
#' @method predict gibclust
#' @exportS3Method
predict.gibclust <- function(object, newdata = NULL, X = NULL, ...) {
  if (is.null(X)) {
    if (is.null(object$training_data)) {
      stop("'X' (the training data) must be supplied. Alternatively, refit with keep_data = TRUE to store training data in the model object.")
    }
    X <- object$training_data
  }
  if (is.null(newdata)) {
    newdata <- X
  }
  if (missing(newdata) || is.null(newdata)) {
    stop("'newdata' must be supplied.")
  }
  if (missing(X) || is.null(X)) {
    stop("'X' (the training data) must be supplied.")
  }
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
  catcols <- object$catcols
  
  if (length(contcols) > 0L && isTRUE(object$scale)) {
    train_means <- colMeans(X[, contcols, drop = FALSE])
    train_sds <- apply(X[, contcols, drop = FALSE], 2, stats::sd)
    X[, contcols] <- scale(X[, contcols, drop = FALSE],
                           center = train_means, scale = train_sds)
    newdata[, contcols] <- scale(newdata[, contcols, drop = FALSE],
                                 center = train_means, scale = train_sds)
  }
  if (length(catcols) > 0L) {
    for (j in catcols) {
      lev <- levels(factor(X[[j]]))
      X[[j]] <- as.integer(factor(X[[j]], levels = lev))
      newdata[[j]] <- as.integer(factor(newdata[[j]], levels = lev))
    }
  }
  
  bws_vec <- numeric(ncol(X))
  if (length(contcols) > 0L) bws_vec[contcols] <- object$s
  if (length(catcols)  > 0L) bws_vec[catcols]  <- object$lambda
  
  eval_list <- coord_to_pxy_eval_R(
    X_train = X,
    X_new = newdata,
    s = if (length(contcols) > 0L) object$s else -1,
    lambda = if (length(catcols)  > 0L) object$lambda else -1,
    cat_cols = catcols,
    cont_cols = contcols,
    contkernel = object$kernels$cont,
    nomkernel = object$kernels$nom,
    ordkernel = object$kernels$ord
  )
  py_x_new <- eval_list$py_x_new
  pxy_train <- coord_to_pxy_R(
    X = X,
    s = if (length(contcols) > 0L) object$s      else -1,
    lambda = if (length(catcols)  > 0L) object$lambda else -1,
    cat_cols = catcols,
    cont_cols = contcols,
    contkernel = object$kernels$cont,
    nomkernel = object$kernels$nom,
    ordkernel = object$kernels$ord
  )
  py_x_train <- pxy_train$py_x
  px_train <- pxy_train$px
  
  qy_t <- .reconstruct_qy_t(py_x_train, object$Cluster, px_train, object$ncl)
  n_new <- ncol(py_x_new)
  ncl <- object$ncl
  kl_mat <- matrix(0, nrow = ncl, ncol = n_new)
  log_p_new <- log(py_x_new)
  log_p_new[!is.finite(log_p_new)] <- 0
  
  for (t in seq_len(ncl)) {
    log_qy_t <- log(qy_t[, t])
    log_qy_t[!is.finite(log_qy_t)] <- 0
    diff_log <- log_p_new - log_qy_t
    kl_mat[t, ] <- colSums(py_x_new * diff_log)
  }
  if (isTRUE(all.equal(object$alpha, 0))) {
    return(apply(kl_mat, 2, which.min))
  }
  qt <- as.numeric(table(factor(object$Cluster, levels = seq_len(ncl))) / object$n)
  if (is.matrix(object$Cluster)) {
    qt <- as.numeric(object$Cluster %*% as.numeric(px_train))
  }
  beta_val <- if (length(object$beta) > 1L) tail(object$beta, 1) else object$beta
  log_qt <- log(qt)
  log_qt[!is.finite(log_qt)] <- -Inf
  log_qt_x <- sweep(-beta_val * kl_mat, 1, log_qt, "+")
  log_qt_x <- sweep(log_qt_x, 2, apply(log_qt_x, 2, max), "-")
  qt_x_new <- exp(log_qt_x)
  qt_x_new <- sweep(qt_x_new, 2, colSums(qt_x_new), "/")
  
  qt_x_new
}