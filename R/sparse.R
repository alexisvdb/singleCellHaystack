#' Returns a row of a sparse matrix of class lgRMatrix.
#' Function made by Ben Bolker and Ott Toomet (see https://stackoverflow.com/questions/47997184/)
#'
#' @param m a sparse matrix of class lgRMatrix
#' @param i the index of the row to return
#'
#' @return A row (logical vector) of the sparse matrix
extract_row_lgRMatrix <- function(m, i = 1) {
  r <- logical(ncol(m)) ## set up vector with FALSE values for results
  inds <- seq(
    from = m@p[i] + 1,
    to = m@p[i + 1],
    length.out = max(0, m@p[i + 1] - m@p[i])
  )
  r[m@j[inds] + 1] <- m@x[inds] ## set values
  return(r)
}

#' Returns a row of a sparse matrix of class dgRMatrix.
#' Function made by Ben Bolker and Ott Toomet (see https://stackoverflow.com/questions/47997184/)
#'
#' @param m a sparse matrix of class dgRMatrix
#' @param i the index of the row to return
#'
#' @return A row (numerical vector) of the sparse matrix
extract_row_dgRMatrix <- function(m, i = 1) {
  r <- numeric(ncol(m)) ## set up vector with zero values for results
  inds <- seq(
    from = m@p[i] + 1,
    to = m@p[i + 1],
    length.out = max(0, m@p[i + 1] - m@p[i])
  )
  r[m@j[inds] + 1] <- m@x[inds] ## set values
  return(r)
}
