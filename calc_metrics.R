calc_metrics <- function(beta, qt, qy_t, hy, quiet = TRUE){
  # Calculates IB performance metrics.
  ht <- entropy(qt)
  hy_t <- crossprod(qt, entropy(qy_t))
  iyt <- hy - hy_t
  L <- ht - beta * iyt
  if (!quiet){
    cat('H(T) =', ht, '\t H(Y|T) =', hy_t, '\t I(Y,T) =', iyt, '\t L =', L, '\n')
  }
  return(list(L, ht, iyt))
}
