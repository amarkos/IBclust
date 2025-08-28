#' Methods for gibclust objects
#'
#' S3 methods available for \code{"gibclust"} objects: \code{print()}, \code{summary()},
#' and \code{plot()}.
#'
#' @name gibclust-methods
#' @aliases print.gibclust summary.gibclust print.summary.gibclust plot.gibclust
#' @keywords methods
#' @seealso \code{\link{DIBmix}}, \code{\link{IBmix}}, \code{\link{GIBmix}}
#' @importFrom graphics barplot points
NULL

# ---- S3: print() -------------------------------------------------------------

#' @rdname gibclust-methods
#' @export
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
  
  mi   <- x$MutualInfo
  ent  <- x$Entropy
  cent <- x$CondEntropy
  
  cat(sprintf("Mutual information I(Y;T): %s\n",
              if (is.finite(mi)) sprintf("%.6f", mi) else "Inf"))
  cat(sprintf("Entropy H(T): %s   Conditional entropy H(T|X): %s\n",
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

# ---- S3: summary() -----------------------------------------------------------

#' @rdname gibclust-methods
#' @export
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
#' @export
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

# ---- S3: plot() --------------------------------------------------------------
#' Plot a gibclust object
#'
#' @param x A gibclust object.
#' @param type Plot type: "sizes" (cluster sizes), "info" (information metrics),
#'   or "beta" (log(beta) trajectory; DIBmix only).
#' @param main Optional title.
#' @param ... Additional arguments passed to base plotting functions.
#' @rdname gibclust-methods
#' @export
plot.gibclust <- function(x, type = c("sizes", "info", "beta"), main = NULL, ...) {
  type <- match.arg(type)
  
  # helper: harden fuzzy memberships for ncl x n membership matrix
  harden <- function(C) apply(C, 2, which.max)
  
  if (type == "sizes") {
    if (isTRUE(all.equal(x$alpha, 0))) {
      # DIBmix: hard labels vector
      cl <- x$Cluster
      if (is.vector(cl) && length(cl) == x$n) {
        tab <- table(as.integer(cl))
        if (is.null(main)) main <- "Cluster sizes (DIBmix)"
        barplot(tab, ylab = "Count", xlab = "Cluster", main = main, ...)
      } else {
        warning("Hard labels unavailable or malformed for DIBmix.")
      }
    } else {
      # IB/GIB: fuzzy membership matrix (ncl x n)
      C <- x$Cluster
      if (is.matrix(C) && nrow(C) == x$ncl && ncol(C) == x$n) {
        lab <- harden(C)
        tab <- table(lab)
        if (is.null(main)) {
          main <- if (isTRUE(all.equal(x$alpha, 1))) "Hardened sizes (IBmix)"
          else "Hardened sizes (GIBmix)"
        }
        barplot(tab, ylab = "Count", xlab = "Cluster", main = main, ...)
      } else {
        warning("Membership matrix unavailable or malformed for IB/GIB.")
      }
    }
    
  } else if (type == "info") {
    vals <- c(`H(T)` = x$Entropy,
              `H(T|X)` = x$CondEntropy,
              `I(Y;T)`    = x$MutualInfo)
    vals[!is.finite(vals)] <- NA_real_
    if (is.null(main)) main <- "Information summary"
    barplot(vals, ylab = "Value", main = main, ...)
    
  } else { # type == "beta"
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
    
    if (is.null(main)) main <- expression(log(beta) ~ " trajectory (DIBmix)")
    
    if (all(is.finite(logb))) {
      plot(idx, logb, type = "l",
           xlab = "Iteration", ylab = expression(log(beta)),
           main = main, ...)
      points(idx, logb, ...)
    } else {
      warning("Non-finite values in log(beta); some points omitted.")
      plot(idx, logb, type = "l",
           xlab = "Iteration", ylab = expression(log(beta)),
           main = main, ...)
      points(idx, logb, ...)
    }
  }
  invisible(x)
}