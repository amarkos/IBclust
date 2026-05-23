# Helper function to convert aibclust merge encoding to hclust convention
.aib_to_hclust_merge <- function(merges, n) {
  m <- nrow(merges)
  hc_merge <- matrix(0L, m, 2)
  slot_event <- -seq_len(n)
  for (step in seq_len(m)) {
    i <- merges[step, 1]
    j <- merges[step, 2]
    hc_merge[step, 1] <- slot_event[i]
    hc_merge[step, 2] <- slot_event[j]
    slot_event[i] <- step
  }
  hc_merge
}

# Helper function to compute leaf order via tree traversal
.aib_compute_order <- function(hc_merge, n) {
  m <- nrow(hc_merge)
  leaves_at_step <- vector("list", m)
  for (step in seq_len(m)) {
    left <- hc_merge[step, 1]
    right <- hc_merge[step, 2]
    left_leaves <- if (left  < 0) -left  else leaves_at_step[[left]]
    right_leaves <- if (right < 0) -right else leaves_at_step[[right]]
    leaves_at_step[[step]] <- c(left_leaves, right_leaves)
  }
  leaves_at_step[[m]]
}

#' Convert an aibclust object to hclust
#'
#' Enables use of standard hierarchical-clustering methods
#' (e.g., \code{\link[stats]{cutree}}, \code{\link[stats]{as.dendrogram}})
#' on AIBmix output.
#'
#' @param x An \code{aibclust} object.
#' @param ... Ignored.
#'
#' @return An object of class \code{"hclust"}.
#'
#' @method as.hclust aibclust
#' @exportS3Method
as.hclust.aibclust <- function(x, ...) {
  n <- x$n
  hc_merge <- .aib_to_hclust_merge(x$merges, n)
  height <- cummax(x$merge_costs)
  ord <- .aib_compute_order(hc_merge, n)
  structure(
    list(
      merge = hc_merge,
      height = height,
      order = ord,
      labels = x$obs_names,
      method = "aib",
      call = match.call(),
      dist.method = "information-bottleneck"
    ),
    class = "hclust"
  )
}

#' @method as.dendrogram aibclust
#' @exportS3Method
as.dendrogram.aibclust <- function(object, ...) {
  as.dendrogram(as.hclust(object), ...)
}