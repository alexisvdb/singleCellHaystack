#' The main Haystack function.
#'
#' @param x a matrix or other object from which coordinates of cells can be extracted.
#' @param dim1 column index or name of matrix for x-axis coordinates.
#' @param dim2 column index or name of matrix for y-axis coordinates.
#' @param assay name of assay data for Seurat method.
#' @param slot name of slot for assay data for Seurat method.
#' @param coord name of coordinates slot for specific methods.
#' @param cutoff cutoff for detection.
#' @param method choose between highD (default) and 2D haystack.
#' @param ... further paramters passed to haystack_raw().
#'
#' @return An object of class "haystack"
#' @export
#'
haystack <- function(x, ...) {
  UseMethod("haystack")
}

#' @rdname haystack
#' @export
haystack.matrix <- function(x, dim1 = 1, dim2 = 2, method = "highD", ...) {
  if (method == "highD")
    haystack_highD(x, ...)
  else
    haystack_2D(x[, dim1], x[, dim2], ...)
}

#' @rdname haystack
#' @export
haystack.data.frame <- function(x, ...) {
  haystack(as.matrix(x), ...)
}

#' @rdname haystack
#' @export
haystack.Seurat <- function(x, assay = "RNA", slot = "data", coord = "tsne", cutoff = 1, ...) {
  y <- GetAssayData(x, slot = slot, assay = assay)
  z <- Embeddings(x, coord)
  haystack(as.matrix(z), detection = y > cutoff, ...)
}

#' @rdname haystack
#' @export
haystack.SingleCellExperiment <- function(x, assay = "counts", coord = "TSNE", cutoff = 1, ...) {
  y <- assay(x, assay)
  z <- reducedDim(x, coord)
  haystack(as.matrix(z), detection = y > cutoff, ...)
}

#' plot_gene_haystack
#'
#' @param x a matrix or other object from which coordinates of cells can be extracted.
#' @param dim1 column index or name of matrix for x-axis coordinates.
#' @param dim2 column index or name of matrix for y-axis coordinates.
#' @param assay name of assay data for Seurat method.
#' @param slot name of slot for assay data for Seurat method.
#' @param coord name of coordinates slot for specific methods.
#' @param ... further paramters passed to plot_gene_haystack_raw().
#'
#' @export
#'
plot_gene_haystack <- function(x, ...) {
  UseMethod("plot_gene_haystack")
}

#' @rdname plot_gene_haystack
#' @export
plot_gene_haystack.matrix <- function(x, dim1 = 1, dim2 = 2, ...) {
  plot_gene_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' @rdname plot_gene_haystack
#' @export
plot_gene_haystack.data.frame <- function(x, dim1 = 1, dim2 = 2, ...) {
  plot_gene_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' @rdname plot_gene_haystack
#' @export
plot_gene_haystack.SingleCellExperiment <- function(x, dim1 = 1, dim2 = 2, assay = "counts", coord = "TSNE", ...) {
  y <- assay(x, assay)
  z <- reducedDim(x, coord)
  plot_gene_haystack_raw(z[, dim1], z[, dim2], expression = y, ...)
}

#' @rdname plot_gene_haystack
#' @export
plot_gene_haystack.Seurat <- function(x, dim1 = 1, dim2 = 2, assay = "RNA", slot = "data", coord = "tsne", ...) {
  y <- GetAssayData(x, slot = slot, assay = assay)
  z <- Embeddings(x, coord)
  plot_gene_haystack_raw(z[, dim1], z[, dim2], expression = y, ...)
}

#' plot_gene_set_haystack
#'
#' @param x a matrix or other object from which coordinates of cells can be extracted.
#' @param dim1 column index or name of matrix for x-axis coordinates.
#' @param dim2 column index or name of matrix for y-axis coordinates.
#' @param assay name of assay data for Seurat method.
#' @param slot name of slot for assay data for Seurat method.
#' @param coord name of coordinates slot for specific methods.
#' @param ... further paramters passed to plot_gene_haystack_raw().
#'
#' @export
#'
plot_gene_set_haystack <- function(x, ...) {
  UseMethod("plot_gene_set_haystack")
}

#' @rdname plot_gene_set_haystack
#' @export
plot_gene_set_haystack.matrix <- function(x, dim1 = 1, dim2 = 2, ...) {
  plot_gene_set_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' @rdname plot_gene_set_haystack
#' @export
plot_gene_set_haystack.data.frame <- function(x, dim1 = 1, dim2 = 2, ...) {
  plot_gene_set_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' @rdname plot_gene_set_haystack
#' @export
plot_gene_set_haystack.SingleCellExperiment <- function(x, dim1 = 1, dim2 = 2, assay = "counts", coord = "TSNE", ...) {
  y <- assay(x, assay)
  z <- reducedDim(x, coord)
  plot_gene_set_haystack_raw(z[, dim1], z[, dim2], detection = y > 1, ...)
}

#' @rdname plot_gene_set_haystack
#' @export
plot_gene_set_haystack.Seurat <- function(x, dim1 = 1, dim2 = 2, assay = "RNA", slot = "data", coord = "tsne", ...) {
  y <- GetAssayData(x, slot = slot, assay = assay)
  z <- Embeddings(x, coord)
  plot_gene_set_haystack_raw(z[, dim1], z[, dim2], detection = y > 1, ...)
}


#' hclust_haystack
#'
#' @param x a matrix or other object from which coordinates of cells can be extracted.
#' @param dim1 column index or name of matrix for x-axis coordinates.
#' @param dim2 column index or name of matrix for y-axis coordinates.
#' @param ... further paramters passed down to methods.
#'
#' @export
#'
hclust_haystack <- function(x, ...) {
  UseMethod("hclust_haystack")
}

#' @rdname hclust_haystack
#' @export
hclut_haystack.matrix <- function(x, dim1 = 1, dim2 = 2, ...) {
  hclust_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' @rdname hclust_haystack
#' @export
hclut_haystack.data.frame <- function(x, dim1 = 1, dim2 = 2, ...) {
  hclust_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' kmeans_haystack
#'
#' @param x a matrix or other object from which coordinates of cells can be extracted.
#' @param dim1 column index or name of matrix for x-axis coordinates.
#' @param dim2 column index or name of matrix for y-axis coordinates.
#' @param ... further paramters passed down to methods.
#'
#' @export
#'
kmeans_haystack <- function(x, ...) {
  UseMethod("kmeans_haystack")
}

#' @rdname kmeans_haystack
#' @export
kmeans_haystack.matrix <- function(x, dim1 = 1, dim2 = 2, ...) {
  kmeans_haystack_raw(x[, dim1], x[, dim2], ...)
}

#' @rdname hclust_haystack
#' @export
kmeans_haystack.data.frame <- function(x, dim1 = 1, dim2 = 2, ...) {
  kmeans_haystack_raw(x[, dim1], x[, dim2], ...)
}
