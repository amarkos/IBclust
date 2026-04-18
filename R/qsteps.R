qy_t_step <- function(py_x, qt_x, qt, px) {
  if (inherits(py_x, "nystrom_matrix")) {
    # Use Nyström factored version
    qy_t <- qy_t_step_nystrom_cpp(
      B = py_x$B,
      col_sums = py_x$col_sums,
      qt_x = qt_x,
      qt = qt,
      px = px
    )
  } else {
    # Use regular dense matrix version
    qy_t <- qy_t_step_cpp(py_x, qt_x, qt, px)
  }
  return(qy_t)
}

qt_x_step_beta <- function(n_rows, T, py_x, qy_t, qt, qt_x) {
  if (inherits(py_x, "nystrom_matrix")) {
    # Use Nyström version
    return(qt_x_step_beta_nystrom_cpp(
      n_rows = n_rows,
      T = T,
      B = py_x$B,
      col_sums = py_x$col_sums,
      qy_t = qy_t,
      qt = qt,
      qt_x = qt_x
    ))
  } else {
    # Use regular version
    return(qt_x_step_beta_cpp(
      n_rows = n_rows,
      T = T,
      py_x = py_x,
      qy_t = qy_t,
      qt = qt,
      qt_x = qt_x
    ))
  }
}

qt_x_step_gib <- function(n_rows, T, beta, alpha, py_x, qy_t, qt) {
  if (inherits(py_x, "nystrom_matrix")) {
    # Use Nystrom version
    return(qt_x_step_gib_nystrom_cpp(
      n_rows = n_rows,
      T = T,
      beta = beta,
      alpha = alpha,
      B = py_x$B,
      col_sums = py_x$col_sums,
      qy_t = qy_t,
      qt = qt
    ))
  } else {
    # Use regular version
    return(qt_x_step_gib_cpp(
      n_rows = n_rows,
      T = T,
      beta = beta,
      alpha = alpha,
      py_x = py_x,
      qy_t = qy_t,
      qt = qt
    ))
  }
}