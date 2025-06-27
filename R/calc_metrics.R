calc_metrics <- function(beta, qt, qy_t, hy, px, qt_x, quiet = TRUE){
  # Calculates IB performance metrics.
  ht <- entropy(qt)
  hy_t <- crossprod(qt, entropy(qy_t))
  iyt <- hy - hy_t
  L_DIB <- ht - beta * iyt
  ht_x <- crossprod(px, entropy(qt_x))
  ixt <- ht - ht_x
  L_IB <- ixt - beta * iyt
  if (!quiet){
    message('H(T) = ', ht, ', H(Y|T) = ', hy_t, ', I(Y,T) = ', iyt,
            ', I(X,T) = ', ixt, ', L_DIB = ', L_DIB, ', L_IB = ', L_IB, '\n')
  }
  return(list(L_DIB, ht, iyt, ixt, L_IB))
}
