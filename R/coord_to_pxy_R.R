#' Compute pxy and Entropies
#'
#' This function computes the pxy matrix and associated entropies using kernel density estimates.
#'
#' @param X Data frame containing the dataset.
#' @param s Bandwidth parameter for continuous variables.
#' @param cat_cols Indices of categorical columns in X.
#' @param cont_cols Indices of continuous columns in X.
#' @param lambda Bandwidth parameter for categorical variables.
#' @param contkernel Continuous kernel (Gaussian or Epanechnikov)
#' @param nomkernel Unordered categorical (nominal) kernel (Aitchison & Aitken or Li & Racine)
#' @param ordkernel Ordered categorical (ordinal) kernel (Li & Racine or Wang & van Ryzin)
#'
#' @return A list containing py_x, px, pxy, hy, hx, hy_x, and mutual_information.
#'
#' @export

# Source the C++ code
#sourceCpp("src/entropy_functions.cpp")

coord_to_pxy_R <- function(X, s, cat_cols, cont_cols, lambda,
                           contkernel = "gaussian",
                           nomkernel = "aitchisonaitken",
                           ordkernel = "liracine"){
  bws <- rep(NA, ncol(X))
  bws[cat_cols] <- lambda
  bws[cont_cols] <- s

  py_x <- t(np::npksum(bws = bws,
                       txdat = X,
                       exdat = X,
                       ckertype = contkernel,
                       ukertype = nomkernel,
                       okertype = ordkernel,
                       return.kernel.weights=TRUE)$kw)
  # Remove grid points with zero density
  if (length(which(rowSums(py_x)==0)) > 0){
    py_x <- py_x[-which(rowSums(py_x)==0),]
  }
  py_x <- sweep(py_x, 2, colSums(py_x),'/')
  px <- matrix(1/nrow(X), nrow = nrow(py_x), ncol = nrow(X))
  pxy <- t(py_x * px)
  # cat('Calculating entropies.\n')
  hx <- entropySingle(rowSums(pxy))
  hy <- entropySingle(colSums(pxy))
  hy_x <- rowSums(pxy) %*% entropy(py_x)
  #hy_x <- eigenMapMatMult(rowSums(pxy), entropy(py_x))
  ixy <- hy - hy_x
  px <- matrix(1/nrow(X), nrow = nrow(X), ncol = 1)
  return(list('py_x' = py_x, 'px' = px, 'pxy' = pxy, 'hy' = hy))
}
