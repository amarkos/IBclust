# R/internal_functions.R

#' Internal Function: IBmix_iterate
#'
#' Performs iterative clustering for the IBmix algorithm.
#'
#' @param X A data frame or matrix containing the dataset.
#' @param ncl Number of clusters.
#' @param beta Regularisation parameter beta.
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
IBmix_iterate <- function(X, ncl, beta, randinit,
                          tol, py_x, hy, px, maxiter, bws_vec,
                          contcols, catcols, runs, verbose = FALSE){
  # Source the C++ code
  #  sourceCpp("src/qt_x_step.cpp")
  
  best_clust <- list()
  Loss <- Inf
  best_clust$Cluster <- rep(NA, nrow(X))
  best_clust$InfoXT <- Inf
  best_clust$InfoYT <- Inf
  best_clust$beta <- beta
  best_clust$s <- bws_vec[contcols]
  best_clust$lambda <- bws_vec[catcols]
  best_clust$ixt <- c()
  best_clust$iyt <- c()
  best_clust$losses <- c()
  if (ncl == 1){
    Loss <- 0
    best_clust$Cluster <- rep(1, nrow(X))
    best_clust$InfoXT <- 0
    best_clust$InfoYT <- 0
    best_clust$beta <- beta
    best_clust$ixt <- 0
    best_clust$iyt <- 0
    best_clust$losses <- 0
  } else {
    for (i in c(1:runs)){
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
      qt_x <- qt_x_step_gib_cpp(n_rows = nrow(X), T = qt_list$T, beta = beta, alpha = 1, py_x, qy_t, as.numeric(qt))
      metrics <- calc_metrics(beta = beta, qt, qy_t, hy, px, qt_x, quiet = TRUE)
      Lval <- metrics[[4]] - beta * metrics[[3]]
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
        qt_x <- qt_x_step_gib_cpp(n_rows = nrow(X), T = qt_list$T, beta = beta, alpha = 1, py_x, qy_t, as.numeric(qt))
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
        Lval <- metrics[[4]] - beta * metrics[[3]]
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
      if (Lval < Loss){
        #   best_clust[[1]] <- Lval
        Loss <- Lval
        best_clust[[1]] <- qt_x
        metrics <- calc_metrics(beta = beta, qt, qy_t, hy, px, qt_x, quiet = TRUE)
        best_clust[[2]] <- as.numeric(metrics[[4]])
        best_clust[[3]] <- as.numeric(metrics[[3]])
        best_clust[[4]] <- beta
      }
      metrics <- calc_metrics(beta = beta, qt, qy_t, hy, px, qt_x, quiet = TRUE)
      best_clust$ixt <- c(best_clust$ixt, metrics[[4]])
      best_clust$iyt <- c(best_clust$iyt, metrics[[3]])
      best_clust$losses <- c(best_clust$losses, metrics[[4]] - beta * metrics[[3]])
      if (verbose){
        message('Run ', i, ' complete.\n')
      }
    }
  }
  
  return(best_clust)
}
