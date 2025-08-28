# R/internal_functions.R

#' Internal Function: DIBmix_iterate
#'
#' Performs iterative clustering for the DIBmix algorithm.
#'
#' @param X A data frame or matrix containing the dataset.
#' @param ncl Number of clusters.
#' @param randinit Optional initial cluster assignments.
#' @param tol Tolerance for convergence.
#' @param py_x Conditional probability matrix \( p(y|x) \).
#' @param hy Entropy \( H(Y) \).
#' @param px Probability matrix \( p(x) \).
#' @param maxiter Maximum number of iterations.
#' @param bws_vec Bandwidth vector.
#' @param contcols Indices of continuous columns.
#' @param catcols Indices of categorical columns.
#' @param runs Number of random starts.
#' @param verbose Defaults to FALSE to suppress progress messages. Change to TRUE to print.
#'
#' @return A list containing clustering results.
#'
#' @keywords internal
DIBmix_iterate <- function(X, ncl, randinit,
                           tol, py_x, hy, px, maxiter, bws_vec,
                           contcols, catcols, runs, verbose = FALSE){
  best_clust <- list()
  Loss <- -Inf
  best_clust$Cluster <- rep(NA, nrow(X))
  best_clust$Entropy <- Inf
  best_clust$CondEntropy <- Inf
  best_clust$MutualInfo <- Inf
  best_clust$InfoXT <- Inf
  best_clust$beta <- NA
  best_clust$alpha <- 0
  best_clust$s <- if (length(contcols) == 0) -1 else as.vector(bws_vec[contcols])
  best_clust$lambda <- if (length(catcols) == 0) -1 else as.vector(bws_vec[catcols])
  best_clust$iters <- NA
  best_clust$converged <- NA
  if (ncl == 1){
    Loss <- 0
    best_clust$Cluster <- rep(1, nrow(X))
    best_clust$Entropy <- 0
    best_clust$CondEntropy <- 0
    best_clust$MutualInfo <- 0
    best_clust$InfoXT <- 0
    best_clust$beta <- 1
    best_clust$iters <- 0
    best_clust$converged <- FALSE
  } else {
    pb <- txtProgressBar(style = 3, min = 0, max = runs)
    for (i in c(1:runs)){
      setTxtProgressBar(pb, i)
      beta_vec <- c()
      # Initialize qt_x (randomly)
      qt_x_init <- matrix(0, nrow = ncl, ncol = nrow(X))
      if (is.null(randinit)){
        rand_init <- sample(rep(1:ncl, each = ceiling(nrow(X) / ncl)), size = nrow(X))
      } else {
        rand_init <- randinit
      }
      for (j in 1:ncl) {
        qt_x_init[j, rand_init == j] <- 1
      }
      #####
      qt_list <- qt_step(X, qt_x_init, ptol = tol, quiet =TRUE)
      qt <- qt_list$qt
      qt_x <- qt_list$qt_x
      qy_t <- qy_t_step_cpp(py_x, qt_x, qt, px)
      qt_x_obj <- qt_x_step_beta_cpp(n_rows = nrow(X), T = qt_list$T, py_x, qy_t, as.numeric(qt), qt_x)
      qt_x <- qt_x_obj$qt_x
      beta <- qt_x_obj$beta
      beta_vec <- c(beta_vec, beta)
      metrics <- calc_metrics(beta = beta, qt, qy_t, hy, px, qt_x, quiet = TRUE)
      Lval <- metrics$iyt
      # Initialize variables for convergence checking
      convergence_threshold <- 1e-5  # Set a small threshold for convergence
      max_iterations <- maxiter  # Prevent infinite loops
      iterations <- 0
      change_in_qt_x <- Inf  # Initialize to Inf to ensure the loop starts

      # Run the iterative process with convergence criteria
      while(change_in_qt_x > convergence_threshold && iterations < max_iterations) {
        iterations <- iterations + 1  # Increment iteration counter

        # Store old qt_x for comparison
        old_qt_x <- qt_x

        # Perform the clustering step
        qt_list <- qt_step(X, qt_x, tol, FALSE)
        qt <- qt_list$qt
        qt_x <- qt_list$qt_x
        qy_t <- qy_t_step_cpp(py_x, qt_x, qt, px)
        qt_x_obj <- qt_x_step_beta_cpp(n_rows = nrow(X), T = qt_list$T, py_x, qy_t, as.numeric(qt), qt_x)
        qt_x <- qt_x_obj$qt_x
        beta <- qt_x_obj$beta
        beta_vec <- c(beta_vec, beta)

        if (nrow(qt_x)!=ncl){
          Lval <- -Inf
          change_in_qt_x <- 0
          next
        } else {
          # Calculate metrics or any other necessary step
          change_in_qt_x <- sum(abs(qt_x - old_qt_x))
        }
        metrics <- calc_metrics(beta = beta, qt, qy_t, hy, px, qt_x, quiet = TRUE)
        Lval <- metrics$iyt
      }
      
      converged_run <- (change_in_qt_x <= convergence_threshold)

      if (Lval > Loss){
        Loss <- Lval
        best_clust$Cluster <- apply(qt_x, 2, function(col) which(col == 1))
        metrics <- calc_metrics(beta = beta, qt, qy_t, hy, px, qt_x, quiet = TRUE)
        best_clust$Entropy <- as.numeric(metrics$ht)
        best_clust$CondEntropy <- as.numeric(metrics$ht_x)
        best_clust$MutualInfo <- as.numeric(metrics$iyt)
        best_clust$InfoXT <- as.numeric(metrics$ixt)
        best_clust$beta <- beta_vec
        best_clust$iters <- iterations
        best_clust$converged <- converged_run
      }
      if (verbose){
        message('Run ', i, ' complete.\n')
      }
    }
    close(pb) 
  }

  return(best_clust)
}

txtProgressBar <- function(min = 0, max = 1, initial = 0, char = "=", width = NA, 
                           title, label, style = 1, file = "") 
{
  if (!identical(file, "") && !(inherits(file, "connection") && 
                                isOpen(file))) 
    stop("'file' must be \"\" or an open connection object")
  if (!style %in% 1L:3L) 
    style <- 1
  .val <- initial
  .killed <- FALSE
  .nb <- 0L
  .pc <- -1L
  nw <- nchar(char, "w")
  if (is.na(width)) {
    width <- getOption("width")
    if (style == 3L) 
      width <- width - 10L
    width <- trunc(width/nw)
  }
  if (max <= min) 
    stop("must have 'max' > 'min'")
  up1 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    if (.nb < nb) {
      cat(strrep(char, nb - .nb), file = file)
      flush.console()
    }
    else if (.nb > nb) {
      cat("\r", strrep(" ", .nb * nw), "\r", strrep(char, 
                                                    nb), sep = "", file = file)
      flush.console()
    }
    .nb <<- nb
  }
  up2 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    if (.nb <= nb) {
      cat("\r", strrep(char, nb), sep = "", file = file)
      flush.console()
    }
    else {
      cat("\r", strrep(" ", .nb * nw), "\r", strrep(char, 
                                                    nb), sep = "", file = file)
      flush.console()
    }
    .nb <<- nb
  }
  up3 <- function(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    pc <- round(100 * (value - min)/(max - min))
    if (nb == .nb && pc == .pc) 
      return()
    cat(paste0("\r  |", strrep(" ", nw * width + 6)), file = file)
    cat(paste(c("\r  |", rep.int(char, nb), rep.int(" ", 
                                                    nw * (width - nb)), sprintf("| %3d%%", pc)), collapse = ""), 
        file = file)
    flush.console()
    .nb <<- nb
    .pc <<- pc
  }
  getVal <- function() .val
  kill <- function() if (!.killed) {
    cat("\n", file = file)
    flush.console()
    .killed <<- TRUE
  }
  up <- switch(style, up1, up2, up3)
  up(initial)
  structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}


setTxtProgressBar <- function (pb, value, title = NULL, label = NULL) 
{
  if (!inherits(pb, "txtProgressBar")) 
    stop(gettextf("'pb' is not from class %s", dQuote("txtProgressBar")), 
         domain = NA)
  oldval <- pb$getVal()
  pb$up(value)
  invisible(oldval)
}
