calculate_dist_to_cells <- function(x, y, method = c("ori", "knn"), k = nrow(x)) {
  method <- match.arg(method)
  switch(method,
         "ori" = calculate_dist_to_cells_ori(x, y),
         "knn" = calculate_dist_to_cells_knn(x, y, k)
  )
}

calculate_dist_to_cells_ori <- function(set1, set2) {
  singleCellHaystack:::get_dist_two_sets(set1, set2)
}


calculate_dist_to_cells_knn <- function(x, y, k = nrow(x)) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    warning("Required package FNN not available.")
  }

  knn <- FNN::get.knnx(x, y, k = k)
  m <- matrix(NA_real_, nrow = nrow(x), ncol = nrow(y))

  for (i in seq_len(nrow(y))) {
    m[knn$nn.index[i, ], i] <- knn$nn.dist[i, ]
  }
  m
}
