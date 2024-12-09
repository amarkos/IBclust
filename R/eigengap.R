eigengap <- function(data, contcols, catcols, bw, ncl){
  powerset <- rje::powerSet(c(1:ncol(data)), ncol(data)-2)
  num_lvls_vec <- sapply(data[, catcols, drop = FALSE], function(x) length(unique(x)))
  lambda_max <- (num_lvls_vec-1)/(num_lvls_vec)
  s_max <- 10
  max_bw <- rep(NA, ncol(data))
  max_bw[contcols] <- s_max
  max_bw[catcols] <- lambda_max
  pxy_list <- coord_to_pxy_R(data,
                             s = bw[contcols],
                             cat_cols = catcols,
                             cont_cols = contcols,
                             lambda = bw[catcols])
  py_x <- pxy_list$py_x
  eigens <- eigen(py_x)$values
  eigen_gap <- Re(eigens)[ncl] - Re(eigens)[ncl+1]
  max_eigen <- eigen_gap
  max_inx <- 1
  for (i in 2:length(powerset)){
    bw_new <- bw
    bw_new[powerset[[i]]] <- max_bw[powerset[[i]]]
    s_new <- bw_new[contcols]
    lambda_new <- bw_new[catcols]
    pxy_list <- coord_to_pxy_R(data, s = s_new,
                               cat_cols = catcols,
                               cont_cols = contcols,
                               lambda = lambda_new)
    py_x <- pxy_list$py_x
    eigens <- eigen(py_x)$values
    eigen_gap <- Re(eigens)[ncl] - Re(eigens)[ncl+1]
    #cat(powerset[[i]], '\t:', eigen_gap, '\n')
    if (eigen_gap > max_eigen){
      max_eigen <- eigen_gap
      max_inx <- i
    }
  }
  if (max_inx == 1){
    return(bw)
  } else {
    bw_new <- bw
    bw_new[powerset[[max_inx]]] <- max_bw[powerset[[max_inx]]]
    return(bw_new)
  }
}
