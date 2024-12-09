entropy <- function(P) {
  if (is.null(dim(P))) { # Vector input
    p <- P
    p_nonzero <- p[p > 0]
    return(-sum(p_nonzero * log2(p_nonzero)))
  } else { # Matrix input
    H <- apply(P, 2, function(col) {
      p_nonzero <- col[col > 0]
      -sum(p_nonzero * log2(p_nonzero))
    })
    return(H)
  }
}
