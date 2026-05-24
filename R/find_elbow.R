#' Detect a knee/elbow in a monotone curve
#'
#' Identifies the point on a curve that lies at maximum perpendicular distance
#' from the straight line connecting two endpoints. For information retention
#' curves this is a common heuristic for choosing the number of clusters.
#'
#' The function expects \code{x} sorted in strictly ascending order, and
#' \code{y} either non-decreasing or non-increasing throughout (no
#' direction reversals). The endpoints \code{x_start} and \code{x_end}
#' must be values present in \code{x}; the segment of the curve between
#' them is searched for the elbow/knee.
#'
#' @param x A numeric vector of x-coordinates, sorted in strictly ascending
#'   order.
#' @param y A numeric vector of y-coordinates the same length as \code{x},
#'   either non-decreasing or non-increasing throughout.
#' @param x_start The x-coordinate of the segment's starting point. Must
#'   appear in \code{x}.
#' @param x_end The x-coordinate of the segment's endpoint. Must appear in
#'   \code{x}, and satisfy \code{x_end > x_start}.
#'
#' @return A list with components:
#'   \item{x}{The x-coordinate of the detected elbow/knee.}
#'   \item{y}{The y-coordinate of the detected elbow/knee.}
#'   \item{index}{The index in the original \code{x}/\code{y} vectors.}
#'   \item{distance}{The perpendicular distance to the connecting line.}
#'
#' @examples
#' # Synthetic monotone-decreasing curve with a clear elbow
#' x <- 1:20
#' y <- 1 / (1 + exp(0.5 * (x - 5)))
#' elbow <- find_elbow(x, y, x_start = 1, x_end = 20)
#' elbow$x
#' plot(x, y, type = "b")
#' points(elbow$x, elbow$y, col = "red", pch = 19, cex = 1.5)
#'
#' # Applied to an AIBmix information retention curve
#' \donttest{
#' fit <- AIBmix(iris[, -5])
#' find_elbow(seq_along(fit$info_ret), fit$info_ret,
#'            x_start = 1, x_end = 10)
#' }
#'
#' @export
find_elbow <- function(x, y, x_start, x_end) {
  # Validation
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("'x' and 'y' must be numeric vectors.")
  }
  if (length(x) != length(y)) {
    stop(sprintf("'x' (length %d) and 'y' (length %d) must have the same length.",
                 length(x), length(y)))
  }
  if (length(x) < 3) {
    stop("Need at least 3 points to detect an elbow.")
  }
  if (is.unsorted(x, strictly = TRUE)) {
    stop("'x' must be sorted in strictly ascending order.")
  }
  dy <- diff(y)
  if (!(all(dy >= 0) || all(dy <= 0))) {
    stop("'y' must be monotone (non-decreasing or non-increasing) throughout.")
  }
  if (missing(x_start) || missing(x_end)) {
    stop("'x_start' and 'x_end' must be supplied.")
  }
  if (length(x_start) != 1 || length(x_end) != 1) {
    stop("'x_start' and 'x_end' must be single numeric values.")
  }
  if (!(x_start %in% x) || !(x_end %in% x)) {
    stop("'x_start' and 'x_end' must both be values present in 'x'.")
  }
  if (x_end <= x_start) {
    stop("'x_end' must be strictly greater than 'x_start'.")
  }
  
  # Subset to the segment between endpoints
  idx_keep <- which(x >= x_start & x <= x_end)
  x_sub <- x[idx_keep]
  y_sub <- y[idx_keep]
  
  if (length(x_sub) < 3) {
    stop("Segment between 'x_start' and 'x_end' must contain at least 3 points.")
  }
  
  # Endpoints of connecting line
  line_start <- c(x_sub[1], y_sub[1])
  line_end <- c(x_sub[length(x_sub)], y_sub[length(y_sub)])
  
  # Perpendicular dist from each interior point to the line
  num <- abs((line_end[2] - line_start[2]) * x_sub -
               (line_end[1] - line_start[1]) * y_sub +
               line_end[1] * line_start[2] -
               line_end[2] * line_start[1])
  denom <- sqrt((line_end[2] - line_start[2])^2 +
                  (line_end[1] - line_start[1])^2)
  distances <- num / denom
  
  elbow_local <- which.max(distances)
  elbow_global <- idx_keep[elbow_local]
  
  list(
    x = x[elbow_global],
    y = y[elbow_global],
    index = elbow_global,
    distance = distances[elbow_local]
  )
}