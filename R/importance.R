# Internal helper: compute I(T; Y_j) for each variable j
#' @keywords internal
#' @noRd
.compute_variable_importance <- function(X, cluster, s, lambda,
                                         contcols, catcols, kernels,
                                         nystrom_landmarks = NULL) {
  n <- nrow(X)
  X <- as.data.frame(X)
  p <- ncol(X)
  
  if (length(cluster) != n) {
    stop("Length of 'cluster' must match nrow(X).")
  }
  cluster <- as.integer(cluster)
  cluster_ids <- sort(unique(cluster))
  ncl <- length(cluster_ids)
  
  # Hard assignment matrix qt_x: ncl x n
  qt_x <- matrix(0, nrow = ncl, ncol = n)
  for (k in seq_len(ncl)) {
    qt_x[k, cluster == cluster_ids[k]] <- 1
  }
  px <- rep(1/n, n)
  qt <- as.vector(qt_x %*% px)
  
  iyt_vec <- numeric(p)
  names(iyt_vec) <- colnames(X)
  
  use_nystrom <- !is.null(nystrom_landmarks)
  
  for (j in seq_len(p)) {
    is_cont <- j %in% contcols
    X_j <- X[, j, drop = FALSE]
    
    if (is_cont) {
      idx <- which(contcols == j)
      s_j <- if (length(s) == 1L) s else s[idx]
      lam_j <- -1
      cat_cols_j <- integer(0)
      cont_cols_j <- 1L
    } else {
      idx <- which(catcols == j)
      s_j <- -1
      lam_j <- if (length(lambda) == 1L) lambda else lambda[idx]
      cat_cols_j <- 1L
      cont_cols_j <- integer(0)
    }
    
    if (use_nystrom) {
      pxy_list <- suppressWarnings(coord_to_pxy_nystrom_R(
        X_j,
        s = s_j,
        lambda = lam_j,
        cat_cols = cat_cols_j,
        cont_cols = cont_cols_j,
        contkernel = kernels$cont,
        nomkernel = kernels$nom,
        ordkernel = kernels$ord,
        n_landmarks = length(nystrom_landmarks),
        landmark_indices = nystrom_landmarks)
      )
    } else {
      pxy_list <- coord_to_pxy_R(
        X_j,
        s = s_j,
        lambda = lam_j,
        cat_cols = cat_cols_j,
        cont_cols = cont_cols_j,
        contkernel = kernels$cont,
        nomkernel = kernels$nom,
        ordkernel = kernels$ord
      )
    }
    
    py_x_j <- pxy_list$py_x
    hy_j <- pxy_list$hy
    
    qy_t_j <- qy_t_step(py_x_j, qt_x, qt, px)
    hy_t_j <- as.numeric(crossprod(qt, entropy(qy_t_j)))
    iyt_vec[j] <- hy_j - hy_t_j
  }
  
  iyt_vec
}


# Internal helper: plot variable importance using base graphics
#' @keywords internal
#' @importFrom graphics par barplot legend
#' @noRd
.plot_variable_importance <- function(iyt, X, color_by_type = TRUE,
                                      col = NULL, main = NULL, ...) {
  # Ascending sort so largest ends up at top in horizontal barplot
  iyt_sorted <- iyt[order(iyt, decreasing = FALSE)]
  
  type_colors <- c(Continuous = "blue",
                   Nominal = "red",
                   Ordinal = "green4")
  
  if (!is.null(col)) {
    cols <- rep_len(col, length(iyt_sorted))
    show_legend <- FALSE
  } else if (color_by_type) {
    var_types <- vapply(X[names(iyt_sorted)], function(column) {
      if (is.ordered(column))      "Ordinal"
      else if (is.factor(column))  "Nominal"
      else                         "Continuous"
    }, character(1))
    cols <- unname(type_colors[var_types])
    show_legend <- TRUE
  } else {
    cols <- "gray40"
    show_legend <- FALSE
  }
  
  if (is.null(main)) main <- "Variable Importance"
  
  # Left margin sized for variable name lengths
  name_width <- max(nchar(names(iyt_sorted)))
  old_par <- graphics::par(
    mar = c(5, min(20, max(6, name_width * 0.45)), 4, 2) + 0.1
  )
  on.exit(graphics::par(old_par), add = TRUE)
  
  graphics::barplot(
    iyt_sorted,
    horiz = TRUE,
    las = 1,
    col = cols,
    border = "black",
    xlab = expression("Mutual Information " * I(T * ";" ~ Y[j])),
    main = main,
    ...
  )
  
  if (show_legend) {
    used_types <- unique(var_types)
    graphics::legend("bottomright",
                     legend = used_types,
                     fill = unname(type_colors[used_types]),
                     bty = "n",
                     title = "Variable Type")
  }
  
  invisible(iyt_sorted)
}