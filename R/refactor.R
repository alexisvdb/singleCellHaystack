# Refactor haystack_continuous_highD.
haystack_continuous_v2 <- function(x, expression, grid.points = 100, weights.advanced.Q = NULL,
                                   scale = TRUE, grid.method = "centroid",
                                   randomization.count = 100, n.genes.to.randomize = 100,
                                   selection.method.genes.to.randomize = "heavytails", grid.coord = NULL,
                                   spline.method = "ns", verbose = TRUE) {
  renderPB <- isFALSE(getOption("rstudio.notebook.executing"))

  message("### calling haystack_continuous_highD()...")

  if (!is.null(grid.coord)) {
    grid.points <- nrow(grid.coord)

    if (ncol(x) != ncol(grid.coord)) {
      stop("Coordinates and grid points have different number of columns.")
    }

    if (!is.null(colnames(x)) || !is.null(colnames(grid.coord))) {
      if (!identical(colnames(x), colnames(grid.coord))) {
        stop("Coordinates and grid points have different column names.")
      }
    }
  }

  # Check for sparseMatrixStats package.
  if (requireNamespace("sparseMatrixStats", quietly = TRUE)) {
    message("### Using package sparseMatrixStats to speed up statistics in sparse matrices.")
    useSMS <- TRUE
  } else {
    message("### Package sparseMatrixStats not found. Install for speed improvements.")
    useSMS <- FALSE
  }

  message("### Calculating row-wise mean, stdev and CV ... ")

  exprs_cv <- calculate_cv(expression)
  sel_bad <- exprs_cv$sd == 0

  expression <- expression[!sel_bad, ]
  exprs_cv <- exprs_cv[!sel_bad, ]

  if (sum(sel_bad) > 0) {
    message("### Filtered ", sum(sel_bad), " genes with zero variance...")
  }

  message("### Using ", randomization.count, " randomizations...")
  message("### Using ", n.genes.to.randomize, " genes to randomize...")

  # check input
  if (!is.numeric(x) && !all(apply(x, 2, is.numeric))) {
    stop("'x' must be a numeric matrix or data.frame")
  }
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a numeric matrix or data.frame")
  }
  # if(ncol(x) < 2)
  #  stop("'x' must have at least 2 columns")
  # if(!is.matrix(expression) && ! inherits(expression, "dgCMatrix") && ! inherits(expression, "dgRMatrix"))
  #    stop("'expression' must be a matrix, dgCMatrix, or dgRMatrix")
  if (ncol(expression) != nrow(x)) {
    stop("The number of columns in 'expression' must be the same as the rows in 'x'")
  }
  if (!is.numeric(grid.points)) {
    stop("The value of 'grid.points' must be a numeric")
  }
  if (grid.points >= ncol(expression)) {
    stop("The number of grid points appears to be very high (higher than the number of cells). You can set the number of grid points using the 'grid.points' parameter.")
  }
  if (!is.null(weights.advanced.Q)) {
    if (!is.numeric(weights.advanced.Q)) {
      stop("'weights.advanced.Q' must either be NULL or a numeric vector")
    }
    if (length(weights.advanced.Q) != nrow(x)) {
      stop("The length of 'weights.advanced.Q' must be the same as the number of rows in 'x'")
    }
  }
  if (spline.method != "ns" & spline.method != "bs") {
    stop("'spline.method' should be either 'ns' or 'bs'")
  }

  count.cells <- ncol(expression)
  count.genes <- nrow(expression)
  if (n.genes.to.randomize > count.genes) {
    warning("Number of genes to randomize (", n.genes.to.randomize, ") is higher than the number of genes in the data.")
    warning("Setting number of genes to randomize to ", count.genes)
    n.genes.to.randomize <- count.genes
  }

  # warn about unusal input sizes
  if (nrow(x) < 50) {
    warning("The number of cells seems very low (", nrow(x), "). Check your input.")
  }
  if (nrow(expression) < 100) {
    warning("The number of genes seems very low (", nrow(expression), "). Check your input.")
  }

  # warn about extreme values for 'grid.points'
  if (grid.points < 10) {
    warning("The value of 'grid.points' appears to be very low (<10). You can set the number of grid points using the 'grid.points' parameter.")
  }
  if (grid.points > count.cells / 10) {
    warning("The value of 'grid.points' appears to be very high (> No. of cells / 10). You can set the number of grid points using the 'grid.points' parameter.")
  }

  # scale data if needed
  if (scale) {
    message("### scaling input data...")
    x <- scale(x)
    # save the mean and stdev of the scaling
    x.scale.center <- attr(x = x, which = "scaled:center")
    x.scale.scale <- attr(x = x, which = "scaled:scale")
  }

  ### Calculate reference distribution Q
  # 1. Using all points, calculate grid.points filling the space.
  # 2. Calculate distance of grid points to all points, and density.
  # 3. Use density to estimate Q.
  # 4. Normalize to sum to 1.
  if (is.null(grid.coord)) {
    message("### deciding grid points...")
    grid.coord <- calculate_grid_points(x, method = grid.method, ngrid_points = grid.points)

    # add another warning for the case that the number of grid.points was changed
    if (nrow(grid.coord) != grid.points) {
      warning("The number of grid points was changed from ", grid.points, " to ", nrow(grid.coord))
      grid.points <- nrow(grid.coord)
    }
  }

  dist.to.grid <- calculate_dist_to_cells(x, grid.coord)

  density.contributions <- calculate_density(dist.to.grid)

  Q <- calculate_Q_dist(density.contributions, weights = weights.advanced.Q, verbose = verbose)

  # Calculate P distribution.
  # 1. Calculate for each feature.
  # 2. Use the density of points around grid.points weighted by the value of each
  # feature.
  P <- calculate_P_dist(expression, density.contributions)

  # Calculate KLD.
  # 1. Use P and Q distributions.
  KLD <- calculate_KLD(P, Q)

  message("### performing randomizations...")
  genes_to_randomize <- select_genes_to_randomize(exprs_cv, n = n.genes.to.randomize, method = selection.method.genes.to.randomize)

  KLD_rand <- randomize_KLD(density.contributions, expression[genes_to_randomize$result$index, ], Q)


  message("### estimating p-values...")
  cv <- exprs_cv$cv
  names(cv) <- rownames(exprs_cv)

  p.vals <- calculate_Pval(KLD, KLD_rand, cv, cv[genes_to_randomize$result$gene])

  rand_info <- p.vals$info
  rand_info$method <- spline.method
  rand_info$genes_to_randomize <- genes_to_randomize
  rand_info$KLD_rand <- KLD_rand
  p.vals <- p.vals$fitted

  # bonferroni correction for multiple testing
  p.adjs <- p.vals + log10(length(p.vals))
  p.adjs[p.adjs > 0] <- 0 # p values should be at most 1; so log10 should be <= 0

  # prepare grid coordinates to return
  # if input data was scaled, the grid points have to be re-scaled
  # else nothing has to be done
  if (scale) {
    grid.coord <- grid.coord * rep(x.scale.scale, each = nrow(grid.coord)) + rep(x.scale.center, each = nrow(grid.coord))
  }

  message("### returning result...")
  # prepare the 'haystack' object to return
  res <- list(
    results = data.frame(
      D_KL = KLD,
      log.p.vals = p.vals,
      log.p.adj = p.adjs,
      row.names = row.names(expression)
    ),
    info = list(
      method = "continuous_v2",
      randomization = rand_info,
      grid.coordinates = grid.coord,
      coord_mean = x.scale.center,
      coord_std = x.scale.scale,
      densities = density.contributions,
      cv = exprs_cv
    )
  )
  class(res) <- "haystack"
  res
}


extract_row <- function(m, i, ...) {
  UseMethod("extract_row")
}

# extract_row.matrix <- function(m, i, drop=TRUE, ...) {
#  m[i, , drop=drop]
# }

extract_row.dgRMatrix <- function(m, i, ...) {
  cut.index1 <- m@p[i]
  cut.index2 <- m@p[i + 1]
  ind <- numeric()
  val <- numeric()
  if (cut.index1 < cut.index2) {
    ind <- m@j[((cut.index1 + 1):cut.index2)] + 1
    val <- m@x[((cut.index1 + 1):cut.index2)]
  }
  list(ind = ind, val = val)
}
# New functions.

calculate_cv <- function(x) {
  message("#> calculate_cv")
  UseMethod("calculate_cv")
}

calculate_cv.matrix <- function(x) {
  mean <- rowMeans(x)
  sd <- apply(x, 1, sd)
  cv <- (sd + 1e-300) / (mean + 1e-300)
  return(data.frame(gene = rownames(x), mean = mean, sd = sd, cv = cv))
}

calculate_cv.dgCMatrix <- function(x) {
  mean <- Matrix::rowMeans(x)
  sd <- sparseMatrixStats::rowSds(x, useNames = FALSE)
  cv <- (sd + 1e-300) / (mean + 1e-300)
  return(data.frame(gene = rownames(x), mean = mean, sd = sd, cv = cv))
}

calculate_cv.dgRMatrix <- function(x) {
  calculate_cv(as(x, "CsparseMatrix"))
}

#' get_grid_points
#'
#' x is a matrix of cell coordinates, either from spatial
#' coordinates or embedding coordinates.
#'
#' ngrid_points specifies the number of grid points.
#'
#' method can be centroid (kmeans) or seeding (kmeans++).
#'
#' Returns matrix of grid_points x cells.
calculate_grid_points <- function(x, method = "centroid", ngrid_points = 100) {
  singleCellHaystack:::get_grid_points(input = x, method = method, grid.points = ngrid_points)
}

#' get_euclidean_distance.
# calculate_euclidean_distance = function(x,y) {
#  singleCellHaystack:::get_euclidean_distance(x, y)
# }

#' calculate_dist_to_cells
#'
#' Calls get_dist_two_sets. Returns the distance from cells to the grid points.
calculate_dist_to_cells <- function(set1, set2) {
  singleCellHaystack:::get_dist_two_sets(set1, set2)
}

#' calculate_density
#'
#' Calculate grid densities.
#'
#' x is a matrix with the distance of cells to the grid
#' points. This returns a matrix with the density of cells around grid points.
#' N x M matrix (N cells x M grid points)
calculate_density <- function(x) {
  bandwidth <- median(apply(x, 1, min))
  x <- x / bandwidth
  exp(-x * x / 2)
}

## Calculate Q.
calculate_Q_dist <- function(x, weights = NULL, verbose = FALSE, pseudo = 1e-300) {
  if (verbose) message("#> Calculating Q distribution ... ", appendLF = FALSE)

  if (is.null(weights)) {
    Q <- colSums(x)
  } else {
    Q <- colSums(x * weights)
  }
  # pseudo <- 1e-300 # quantile(Q[Q>0],0.01)
  Q <- Q + pseudo
  Q <- Q / sum(Q) # normalize to sum to 1

  if (verbose) message("done.", appendLF = TRUE)

  Q
}

#' calculate_P_dist
#'
#' Returns a matrix with the P distribution for each gene x grid point.
#'
#'
calculate_P_dist <- function(expression, density, ...) {
  # message("#> calculate_P_dist")
  UseMethod("calculate_P_dist")
}

calculate_P_dist.matrix <- function(expression, density, pseudo = 1e-300, ...) {
  ngenes <- nrow(expression)
  ngrid_points <- ncol(density)

  P <- matrix(0, nrow = ngenes, ncol = ngrid_points)

  for (i in seq_len(ngenes)) {
    Pi <- colSums(density * expression[i, ])
    Pi <- Pi + pseudo
    Pi <- Pi / sum(Pi)
    P[i, ] <- Pi
  }
  P
}

calculate_P_dist.dgRMatrix <- function(expression, density, pseudo = 1e-300, ...) {
  ngenes <- nrow(expression)
  ngrid_points <- ncol(density)

  P <- matrix(0, nrow = ngenes, ncol = ngrid_points)

  for (i in seq_len(ngenes)) {
    wl <- extract_row(expression, i)
    val <- wl$val
    ind <- wl$ind

    Pi <- colSums(density[ind, , drop = FALSE] * val)
    Pi <- Pi + pseudo
    Pi <- Pi / sum(Pi)
    P[i, ] <- Pi
  }
  P
}

calculate_P_dist.dgCMatrix <- function(expression, density, pseudo = 1e-300, ...) {
  expression <- as(expression, "RsparseMatrix")
  calculate_P_dist(expression, density, pseudo = 1e-300, ...)
}

#' calculate_KLD
#'
#' Takes a matrix of P distributions and the Q distribution and
#' calculates KLD.
calculate_KLD <- function(P, Q) {
  rowSums(P * log(t(t(P) / Q)))
}

#' @param x CV
select_genes_to_randomize <- function(x, n = 100, method = "heavytails", tail = 9, verbose = FALSE) {
  cv <- x$cv
  ngenes <- nrow(x)
  #  coeffVar  <- expr.sd/expr.mean #
  o <- order(cv)

  if (method == "uniform") {
    index <- o[floor(seq(1, nrow(x), length.out = n))]
  } else if (method == "heavytails") {
    # ls <- head(o, tail)
    # rs <- tail(o, tail)
    # ms <- o[seq(tail + 1, nrow(x)-tail-1, length.out=n-tail*2)]
    # index <- c(ls, ms, rs)
    index <- o[c(
      1:9,
      as.integer(seq(10, ngenes - 9, length.out = n - 18)),
      length(o) - (8:0)
    )]
    # index <- o[c(1:9,
    #             as.integer(seq(10, ngenes-9, length.out=n-18)),
    #             length(o)-(8:0))]
  }

  list(data = data.frame(gene = rownames(x), cv = cv), result = data.frame(index = index, gene = rownames(x)[index], cv = cv[index]))
}

#' randomize_KLD
#'
#' Takes as arguments the matrix of densities, the expression matrix
#' of the genes to randomize and the reference distribution Q.
#' Shuffles the columns of the expression matrix in each randomization
#' and calculates KLD with calculate_KLD.
#'
randomize_matrix <- function(x) {
  # message("#> randomize_matrix")
  UseMethod("randomize_matrix")
}

randomize_matrix.matrix <- function(x) {
  x[, sample(ncol(x))] # TODO: support wrswoR::sample_int_expj.
}

randomize_matrix.dgCMatrix <- function(x) {
  x[, sample(ncol(x))] # TODO: support wrswoR::sample_int_expj.
}

randomize_matrix.dgRMatrix <- function(x) {
  x <- as(x, "CsparseMatrix")
  x <- randomize_matrix(x)
  as(x, "RsparseMatrix")
}

randomize_KLD <- function(density, expression, Q, n_randomizations = 100, pseudo = 1e-300, verbose = FALSE) {
  KLD <- matrix(NA_real_,
    nrow = nrow(expression),
    ncol = n_randomizations,
    dimnames = list(rownames(expression), seq_len(n_randomizations))
  )
  for (n in seq_len(n_randomizations)) {
    expression <- randomize_matrix(expression)
    P <- calculate_P_dist(expression, density)
    KLD[, n] <- calculate_KLD(P, Q)
  }
  KLD
}

calculate_Pval <- function(kld, kld_rand, exprs_cv, rand_cv, method = "ns", ...) {
  get_log_p_D_KL_continuous(
    D_KL.observed = kld,
    D_KL.randomized = kld_rand,
    all.coeffVar = exprs_cv,
    train.coeffVar = rand_cv,
    spline.method = method, ...
  )
}
