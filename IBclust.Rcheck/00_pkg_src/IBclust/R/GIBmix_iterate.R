# R/internal_functions.R

#' Internal Function: IBmix_iterate
#'
#' Performs iterative clustering for the IBmix algorithm.
#'
#' @param X A data frame or matrix containing the dataset.
#' @param ncl Number of clusters.
#' @param beta Regularisation parameter beta.
#' @param alpha Strength of relative entropy alpha.
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
GIBmix_iterate <- function(X, ncl, beta, alpha, randinit,
                           tol, py_x, hy, px, maxiter, bws_vec,
                           contcols, catcols, runs, verbose = FALSE){
  # Source the C++ code
  #  sourceCpp("src/qt_x_step.cpp")
  
  best_clust <- list()
  Loss <- Inf
  best_clust$Cluster <- rep(NA, nrow(X))
  best_clust$Entropy <- Inf
  best_clust$CondEntropy <- Inf
  best_clust$MutualInfo <- Inf
  best_clust$InfoXT <- Inf
  best_clust$beta <- beta
  best_clust$alpha <- alpha
  best_clust$s <- bws_vec[contcols]
  best_clust$lambda <- bws_vec[catcols]
  best_clust$ht <- c()
  best_clust$hy_t <- c()
  best_clust$iyt <- c()
  best_clust$losses <- c()
  if (ncl == 1){
    Loss <- 0
    best_clust$Cluster <- rep(1, nrow(X))
    best_clust$Entropy <- 0
    best_clust$CondEntropy <- 0
    best_clust$MutualInfo <- 0
    best_clust$InfoXT <- 0
    best_clust$beta <- beta
    best_clust$alpha <- alpha
    best_clust$ht <- 0
    best_clust$hy_t <- 0
    best_clust$iyt <- 0
    best_clust$losses <- 0
  } else {
    pb <- txtProgressBar(style = 3, min = 0, max = runs)
    for (i in c(1:runs)){
      setTxtProgressBar(pb, i)
      #set.seed(i)
      # 2. Initialize qt_x (randomly)
      qt_x_init <- matrix(0, nrow = ncl, ncol = nrow(X))
      if (is.null(randinit)){
        rand_init <- sample(rep(1:ncl, each = ceiling(nrow(X) / ncl)), size = nrow(X))
      } else {
        rand_init <- randinit
      }
      #if (length(unique(table(rand_init))) == 1){
      #  level1 <- which(rand_init == 1)[1]
      #  rand_init[level1] <- 2
      #}
      for (j in 1:ncl) {
        qt_x_init[j, rand_init == j] <- 1
      }
      #####
      qt_list <- qt_step(X, qt_x_init, ptol = tol, quiet = TRUE)
      qt <- qt_list$qt
      qt_x <- qt_list$qt_x
      qy_t <- qy_t_step_cpp(py_x, qt_x, qt, px)
      qt_x <- qt_x_step_gib_cpp(n_rows = nrow(X), T = qt_list$T, beta = beta, alpha = alpha, py_x, qy_t, as.numeric(qt))
      metrics <- calc_metrics(beta = beta, qt, qy_t, hy, px, qt_x, quiet = TRUE)
      Lval <- metrics[[1]] - alpha * metrics[[2]] - beta * metrics[[3]]
      #cat('I(Y;T) =', Lval, '\n')
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
        
        # Store old Lval for comparison
        #Lval_old <- Lval
        
        # Perform the clustering step
        qt_list <- qt_step(X, qt_x, tol, FALSE)
        qt <- qt_list$qt
        qt_x <- qt_list$qt_x
        qy_t <- qy_t_step_cpp(py_x, qt_x, qt, px)
        qt_x <- qt_x_step_gib_cpp(n_rows = nrow(X), T = qt_list$T, beta = beta, alpha = alpha, py_x, qy_t, as.numeric(qt))
        #if (sum(qt_x) == 0){
        #  Lval <- -Inf
        #  change_in_qt_x <- 0
        #  message('Bad seed.')
        #  next
        #}
        
        if (nrow(qt_x)!=ncl){
          Lval <- -Inf
          change_in_qt_x <- 0
          next
          #ncl_temp <- nrow(qt_x)
          #change_in_qt_x <- Inf
        } else {
          # Calculate metrics or any other necessary step
          #Lval <- calc_metrics(beta = beta, qt, qy_t, hy, quiet = TRUE)[[1]]
          change_in_qt_x <- sum(abs(qt_x - old_qt_x))
        }
        #Lval <- calc_metrics(beta = beta, qt, qy_t, hy, quiet = TRUE)[[1]]
        metrics <- calc_metrics(beta = beta, qt, qy_t, hy, px, qt_x, quiet = TRUE)
        Lval <- metrics[[1]] - alpha * metrics[[2]] - beta * metrics[[3]]
        ### STOP BASED ON LVAL
        #if (Lval < Lval_old){
        #  qt_x <- old_qt_x
        #  beta_vec <- beta_vec[-length(beta_vec)]
        #  qt_list <- qt_step(X, qt_x, tol, FALSE)
        #  qt <- qt_list$qt
        #  qt_x <- qt_list$qt_x
        #  qy_t <- qy_t_step_cpp(py_x, qt_x, qt, px)
        #  break
        #}
        #cat('I(Y;T) =', Lval, '\n')
        #result_vector <- apply(qt_x, 2, function(col) which(col == 1))
        #return(result_vector)
      }
      
      # Optional: Print the change to monitor progress
      # cat("Iteration:", iterations, "- Change in qt_x:", change_in_qt_x, "\n")
      # Removed conditions: & nrow(qt_x)==ncl & !all(apply(qt_x, 2, function(col) which(col == 1)) == rand_init)
      #if (Lval < best_clust[[1]]){
      if (Lval < Loss & nrow(qt_x)==ncl){
        #   best_clust[[1]] <- Lval
        Loss <- Lval
        best_clust[[1]] <- qt_x
        metrics <- calc_metrics(beta = beta, qt, qy_t, hy, px, qt_x, quiet = TRUE)
        best_clust[[2]] <- as.numeric(metrics[[1]])
        best_clust[[3]] <- as.numeric(metrics[[2]])
        best_clust[[4]] <- as.numeric(metrics[[3]])
        best_clust[[5]] <- as.numeric(metrics[[4]])
        best_clust[[6]] <- beta
        best_clust[[7]] <- alpha
      }
      metrics <- calc_metrics(beta = beta, qt, qy_t, hy, px, qt_x, quiet = TRUE)
      best_clust$ht <- c(best_clust$ht, metrics[[1]])
      best_clust$hy_t <- c(best_clust$hy_t, metrics[[2]])
      best_clust$iyt <- c(best_clust$iyt, metrics[[3]])
      best_clust$losses <- c(best_clust$losses, metrics[[1]] - alpha * metrics[[2]] - beta * metrics[[3]])
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

