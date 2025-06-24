qt_step <- function(X, qt_x, ptol, quiet){
  # Performs q(t) update step for generalized Information Bottleneck.
  px <- rep(1 / nrow(X), nrow(X))
  qt <- as.matrix(qt_x %*% px)
  T <- length(qt)
  dropped <- qt <= ptol  # clusters to drop due to exactly/near-zero prob
  if (any(dropped)) {
    qt <- qt[!dropped] # drop unused clusters
    qt_x <- qt_x[!dropped, ]
    T <- length(qt)  # update number of clusters
    if (T > 1){
      qt_x <- qt_x * matrix(1 / colSums(qt_x), nrow = T, ncol = nrow(X))  # renormalize
    } else {
      qt_x <- qt_x * matrix(1 / qt_x, nrow = T, ncol = nrow(X))
    }
    qt <- as.matrix(qt_x %*% px)

    if (!quiet) {
      message(sprintf('%i cluster(s) dropped. Down to %i cluster(s).\n', sum(dropped), T))
    }
  }

  return(list(qt = qt, qt_x = qt_x, T = T))
}
