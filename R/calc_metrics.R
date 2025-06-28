calc_metrics <- function(beta, qt, qy_t, hy, px, qt_x, quiet = TRUE){
  # Calculates IB performance metrics.
  ht <- entropy(qt)
  hy_t <- crossprod(qt, entropy(qy_t))
  iyt <- hy - hy_t
  ht_x <- crossprod(px, entropy(qt_x))
  ixt <- ht - ht_x
  if (!quiet){
    message('H(T) = ', ht, ', H(Y|T) = ', hy_t, ', I(Y,T) = ', iyt,
            ', I(X,T) = ', ixt, '\n')
  }
  return(list(ht, ht_x, iyt, ixt))
}
