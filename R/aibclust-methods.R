#' Methods for aibclust objects
#'
#' S3 methods available for \code{"aibclust"} objects: \code{print()}, \code{summary()},
#' \code{print.summary()}, and \code{plot()}.
#'
#' @name aibclust-methods
#' @aliases print.aibclust summary.aibclust print.summary.aibclust plot.aibclust
#' @keywords methods
#' @importFrom graphics plot points barplot
NULL

#' @rdname aibclust-methods
#' @export
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
#' @export
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

#' @export
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

#' Plot an aibclust object
#'
#' @param x an \code{aibclust} object.
#' @param type \code{"dendrogram"} (default) or \code{"info"} (information retained curve).
#' @param main optional title.
#' @param labels logical; show labels on dendrogram.
#' @param ... passed to graphics.
#' @rdname aibclust-methods
#' @export
plot.aibclust <- function(x, type = c("dendrogram", "info"), main = NULL, labels = TRUE, ...) {
  type <- match.arg(type)
  
  if (type == "dendrogram") {
    # rebuild dendrogram at plot time
    lab <- if (isTRUE(labels)) x$obs_names else NULL
    d <- make_dendrogram(x$merges, x$merge_costs, labels = lab)
    if (is.null(main)) main <- "AIBmix dendrogram"
    plot(d, main = main, ...)
  } else {  # "info"
    # plot information retained vs number of clusters m
    m <- seq_len(x$n) 
    y <- x$info_ret
    if (is.null(main)) main <- "Information retention curve"
    plot(m, y, type = "l", xlab = "Number of clusters (m)",
         ylab = expression(I(T[m] * ";" * Y) / I(X * ";" * Y)), main = main, ...)
    points(m, y, ...)
  }
  invisible(x)
}