#' Methods for aibclust objects
#'
#' S3 methods available for \code{"aibclust"} objects, including extractors
#' for the cluster assignments and model parameters, an information-metrics
#' accessor, conversion methods to \code{hclust} and \code{dendrogram}, and
#' diagnostic plotting.
#'
#' @details
#' The following methods are available:
#' \itemize{
#'   \item \code{\link[=print.aibclust]{print}} and \code{\link[=summary.aibclust]{summary}}:
#'     concise and detailed descriptions of the cluster hierarchy.
#'   \item \code{\link[=fitted.aibclust]{fitted}}: extract the cluster
#'     partition at a requested number of clusters via the \code{ncl}
#'     argument.
#'   \item \code{\link[=coef.aibclust]{coef}}: extract the model's
#'     bandwidth hyperparameters (\code{s}, \code{lambda}).
#'   \item \code{\link[=info_metrics]{info_metrics}}: extract
#'     information-theoretic quantities. Optional \code{ncl} argument
#'     returns scalar values at the chosen cluster count.
#'   \item \code{\link[=as.hclust.aibclust]{as.hclust}} and
#'     \code{\link[=as.dendrogram.aibclust]{as.dendrogram}}: convert to
#'     standard \code{hclust} and dendrogram objects, enabling
#'     \code{cutree()} and dendrogram-based tools.
#'   \item \code{\link[=plot.aibclust]{plot}}: produce diagnostic plots
#'     (\code{type = "dendrogram"}, \code{"info"}, \code{"importance"}, or \code{"similarity"}).
#' }
#'
#' @name aibclust-methods
#' @aliases print.aibclust summary.aibclust print.summary.aibclust plot.aibclust fitted.aibclust coef.aibclust
#' @seealso \code{\link{AIBmix}}
#' @importFrom graphics plot points barplot
NULL

#' @rdname aibclust-methods
#' @param x an \code{aibclust} object.
#' @param object an \code{aibclust} object.
#' @param ... additional arguments passed to or from other methods.
#' @method print aibclust
#' @exportS3Method
print.aibclust <- function(x, ...) {
  header <- "Hierarchical clustering with AIBmix"
  cat(header, "\n", strrep("-", nchar(header)), "\n", sep = "")
  cat("Call: "); print(x$call)
  cat(sprintf("Observations: %d\n", x$n))
  
  # Merge / info overview
  mc <- x$merge_costs
  cat(sprintf("Merges: %d steps  |  merge cost range: [%.6f, %.6f]\n",
              length(mc), min(mc, na.rm = TRUE), max(mc, na.rm = TRUE)))
  cat(sprintf("I(X;Y): %.6f\n", x$I_X_Y))
  
  # Show info retained at a few representative cluster counts
  k_show <- sort(unique(pmin(x$n, c(2L, 3L, 5L, 10L, x$n))))
  k_show <- k_show[k_show >= 1L & k_show <= x$n]
  info_tbl <- data.frame(m = k_show, info_ret = round(x$info_ret[k_show], 6))
  cat("Information retained I(T_m;Y)/I(X;Y) at selected m:\n")
  print(info_tbl, row.names = FALSE)
  
  # Bandwidths and data makeup
  len_safe <- function(z) if (is.null(z)) 0L else length(z)
  cat(sprintf("Continuous variables: %d   Categorical variables: %d\n",
              len_safe(x$contcols), len_safe(x$catcols)))
  
  .format_vec <- function(v, name, digits = 4, max_show = 6) {
    if (!length(v)) return()
    r <- round(v, digits)
    if (length(r) > max_show) {
      shown <- paste(r[1:max_show], collapse = ", ")
      cat(sprintf("%s = %s, ... (%d total)\n", name, shown, length(r)))
    } else {
      cat(sprintf("%s = %s\n", name, paste(r, collapse = ", ")))
    }
  }
  .format_vec(x$s, "s")
  .format_vec(x$lambda, "lambda")
  
  ks <- x$kernels
  if (is.list(ks)) {
    cat(sprintf("kernels = cont:%s, nom:%s, ord:%s\n", ks$cont, ks$nom, ks$ord))
  }
  invisible(x)
}

#' @rdname aibclust-methods
#' @param k Optional integer vector of cluster counts (cuts) to summarise (same as \code{m}).
#' @param m Optional synonym for \code{k}; the number of clusters.
#' @method summary aibclust
#' @exportS3Method
summary.aibclust <- function(object, k = NULL, m = NULL, ...) {
  n <- object$n
  # vectors (guard against NULL)
  info_ret <- if (!is.null(object$info_ret)) object$info_ret else numeric(0)
  I_T_Y <- if (!is.null(object$I_T_Y)) object$I_T_Y else numeric(0)
  parts <- if (!is.null(object$partitions)) object$partitions else list()
  
  # how many cuts can we actually support?
  max_possible <- min(
    n,
    length(info_ret),
    length(I_T_Y),
    length(parts)
  )
  if (!is.finite(max_possible)) max_possible <- 0L
  
  # choose cuts to show
  if (is.null(k) && is.null(m)) {
    k_show <- sort(unique(pmin(n, c(2L, 3L, 5L, 10L, n))))
  } else {
    k_show <- unique(as.integer(c(k, m)))
    k_show <- k_show[is.finite(k_show)]
  }
  k_show <- k_show[k_show >= 1L & k_show <= n]
  
  # clamp to available range
  if (max_possible >= 1L) {
    k_show <- k_show[k_show <= max_possible]
    if (!length(k_show)) k_show <- max_possible
  } else {
    # nothing available; produce empty summaries
    k_show <- integer(0)
  }
  
  # build info table
  if (length(k_show)) {
    info_by_m <- data.frame(
      m = k_show,
      `I(T[m];Y)/I(X;Y)` = info_ret[k_show],
      `I(T[m];Y)` = I_T_Y[k_show],
      check.names = FALSE
    )
  } else {
    warning("info_ret or I(T;Y) unavailable; producing empty summary.")
    info_by_m <- data.frame(
      m = integer(0),
      `I(T[m];Y)/I(X;Y)` = numeric(0),
      `I(T[m];Y)` = numeric(0),
      check.names = FALSE
    )
  }
  
  # sizes per cut (names must match length(k_show))
  sizes_by_m <- setNames(vector("list", length(k_show)), paste0("m=", k_show))
  if (length(k_show)) {
    for (i in seq_along(k_show)) {
      mm <- k_show[i]
      part <- parts[[mm]]
      if (!is.null(part) && length(part) == n) {
        sizes_by_m[[i]] <- sort(table(as.integer(part)), decreasing = TRUE)
      } else {
        sizes_by_m[[i]] <- integer(0)
      }
    }
  }
  
  out <- list(
    call = object$call,
    n = n,
    merges = if (n > 0) n - 1L else 0L,
    I_X_Y = if (length(object$I_X_Y)) object$I_X_Y else NA_real_,
    merge_cost_summary = if (length(object$merge_costs)) {
      c(
        min = min(object$merge_costs, na.rm = TRUE),
        median = stats::median(object$merge_costs, na.rm = TRUE),
        max = max(object$merge_costs, na.rm = TRUE)
      )
    } else c(min = NA_real_, median = NA_real_, max = NA_real_),
    cuts = k_show,
    info_by_m = info_by_m,
    sizes_by_m = sizes_by_m,
    s = object$s,
    lambda = object$lambda,
    contcols = object$contcols,
    catcols = object$catcols,
    kernels = object$kernels
  )
  class(out) <- "summary.aibclust"
  out
}

#' @rdname aibclust-methods
#' @method print summary.aibclust
#' @exportS3Method
print.summary.aibclust <- function(x, ...) {
  header <- "Summary of AIBmix clustering"
  cat(header, "\n", strrep("-", nchar(header)), "\n", sep = "")
  cat("Call: "); print(x$call)
  cat(sprintf("n = %d, merges = %d\n", x$n, x$merges))
  cat(sprintf("I(X;Y) = %.6f\n", x$I_X_Y))
  
  mc <- x$merge_cost_summary
  cat(sprintf("merge cost summary: min=%.6f, median=%.6f, max=%.6f\n",
              mc["min"], mc["median"], mc["max"]))
  
  cat("\nInformation retained by number of clusters m:\n")
  df <- x$info_by_m
  nm_ratio <- "I(T[m];Y)/I(X;Y)"
  nm_abs   <- "I(T[m];Y)"
  
  if (nrow(df)) {
    # make sure they are numeric and round (in case anything became character)
    df[[nm_ratio]] <- round(as.numeric(df[[nm_ratio]]), 6)
    df[[nm_abs  ]] <- round(as.numeric(df[[nm_abs  ]]), 6)
    print(df, row.names = FALSE)
  } else {
    cat("(no information available)\n")
  }
  
  # show a few cluster-size tables
  show_cuts <- head(x$cuts, 3L)
  for (mm in show_cuts) {
    nm <- paste0("m=", mm)
    cat("\nCluster sizes at ", nm, ":\n", sep = "")
    tab <- x$sizes_by_m[[nm]]
    if (length(tab)) print(tab) else cat("(unavailable)\n")
  }
  if (length(x$cuts) > length(show_cuts)) {
    cat(sprintf("\n(Additional cuts summarised: %s)\n",
                paste(paste0("m=", setdiff(x$cuts, show_cuts)), collapse = ", ")))
  }
  
  len_safe <- function(z) if (is.null(z)) 0L else length(z)
  cat(sprintf("\nContinuous variables: %d   Categorical variables: %d\n",
              len_safe(x$contcols), len_safe(x$catcols)))
  
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
  .format_vec(x$s, "s"); .format_vec(x$lambda, "lambda")
  
  ks <- x$kernels
  if (is.list(ks)) cat(sprintf("kernels = cont:%s, nom:%s, ord:%s\n", ks$cont, ks$nom, ks$ord))
  invisible(x)
}

#' @rdname aibclust-methods
#' @param type Plot type: \code{"dendrogram"} (default), \code{"info"} (information retained curve), \code{"importance"} (variable importance bar chart), or \code{"similarity"} (heatmap of the kernel similarity matrix \eqn{P_{Y|X}}).
#' @param main Optional title.
#' @param col Optional color (or, for \code{type = "similarity"}, a colour palette vector).
#' @param labels Logical; show labels on dendrogram.
#' @param X Original data frame used to fit \code{x}; required for
#'   \code{type = "importance"} and \code{type = "similarity"} when the
#'   fit was constructed with \code{keep_data = FALSE}. If the fit already
#'   contains the training data (\code{keep_data = TRUE}), any supplied
#'   \code{X} is ignored with a warning.
#' @param ncl Number of clusters to use. Required for \code{fitted} (cut
#'   the hierarchy) and \code{plot} with \code{type = "importance"}. For
#'   \code{type = "similarity"}, used (if supplied) to reorder rows and
#'   columns of the similarity matrix by the partition at the chosen cut
#'   and to draw cluster-boundary boxes.
#' @param color_by_type Logical; if \code{TRUE}, colour bars by variable type
#'   (continuous / nominal / ordinal). Defaults to \code{TRUE}. Used only by
#'   \code{type = "importance"}.
#' @param colorbar Logical; if \code{TRUE} (default), draw a horizontal
#'   colour scale below the similarity heatmap. Used only by
#'   \code{type = "similarity"}.
#' @method plot aibclust
#' @exportS3Method
plot.aibclust <- function(x, type = c("dendrogram", "info", "importance", "similarity"),
                          X = NULL, ncl = NULL, color_by_type = TRUE, col = NULL,
                          colorbar = TRUE,
                          main = NULL, labels = TRUE, ...) {
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
    if (is.null(ncl)) {
      stop("Argument 'ncl' (number of clusters to cut at) is required for type = 'importance'.")
    }
    if (ncl < 2 || ncl > length(x$partitions)) {
      stop(sprintf("'ncl' must be between 2 and %d.", length(x$partitions)))
    }
    if (nrow(X) != x$n) {
      stop(sprintf("nrow(X) = %d does not match the fitted model's n = %d.",
                   nrow(X), x$n))
    }
    
    cluster <- x$partitions[[ncl]]
    
    iyt <- .compute_variable_importance(
      X = X,
      cluster = cluster,
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
      s = if (length(contcols) > 0L) x$s      else -1,
      lambda = if (length(catcols)  > 0L) x$lambda else -1,
      cat_cols = catcols,
      cont_cols = contcols,
      contkernel = x$kernels$cont,
      nomkernel = x$kernels$nom,
      ordkernel = x$kernels$ord
    )
    M <- pxy_list$py_x
    cluster_sizes <- NULL
    if (!is.null(ncl)) {
      if (ncl < 1 || ncl > length(x$partitions)) {
        stop(sprintf("'ncl' must be between 1 and %d.", length(x$partitions)))
      }
      cl <- as.integer(x$partitions[[ncl]])
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
  
  if (type == "dendrogram") {
    lab <- if (isTRUE(labels)) x$obs_names else NULL
    d <- make_dendrogram(x$merges, x$merge_costs, labels = lab)
    if (is.null(main)) main <- "AIBmix dendrogram"
    if (is.null(col)) {
      plot(d, main = main, ...)
    } else {
      plot(d, main = main, edgePar = list(col = col), ...)
    }
  } else { 
    m <- seq_len(x$n)
    y <- x$info_ret
    line_col <- if (is.null(col)) "black" else col
    if (is.null(main)) main <- "Information retention curve"
    plot(m, y, type = "l", col = line_col,
         xlab = "Number of clusters (m)",
         ylab = expression(I(T[m] * ";" * Y) / I(X * ";" * Y)),
         main = main, ...)
    points(m, y, col = line_col, ...)
  }
  invisible(x)
}

#' @rdname aibclust-methods
#' @method fitted aibclust
#' @exportS3Method
fitted.aibclust <- function(object, ncl, ...) {
  if (missing(ncl) || is.null(ncl)) {
    stop("Number of clusters 'ncl' must be supplied for aibclust objects.")
  }
  if (!is.numeric(ncl) || length(ncl) != 1L || ncl != round(ncl)) {
    stop("'ncl' must be a single integer.")
  }
  ncl <- as.integer(ncl)
  if (ncl < 1L || ncl > length(object$partitions)) {
    stop(sprintf("'ncl' must be between 1 and %d.", length(object$partitions)))
  }
  object$partitions[[ncl]]
}

#' @rdname aibclust-methods
#' @method coef aibclust
#' @exportS3Method
coef.aibclust <- function(object, ...) {
  out <- list()
  if (length(object$contcols) > 0L) out$s <- object$s
  if (length(object$catcols)  > 0L) out$lambda <- object$lambda
  out
}