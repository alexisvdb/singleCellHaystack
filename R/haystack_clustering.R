

### hclust_haystack
#' Function for hierarchical clustering of genes according to their distribution on a 2D plot.
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param detection A logical matrix showing which genes (rows) are detected in which cells (columns)
#' @param genes A set of genes (of the 'detection' data) which will be clustered.
#' @param method The method to use for hierarchical clustering. See '?hclust' for more information. Default: "ward.D".
#'
#' @return An object of class hclust, describing a hierarchical clustering tree.
#' @export
#'
#' @examples
#' # using the toy example of the singleCellHaystack package
#' # define a logical matrix with detection of each gene (rows) in each cell (columns)
#' dat.detection <- dat.expression > 1
#'
#' # running haystack in default mode
#' res <- haystack(dat.tsne, detection=dat.detection, method = "2D")
#'
#' # get biased genes, store in variable gene.subset
#' sorted.table <- show_result_haystack(res, p.value.threshold = 1e-5)
#' gene.subset <- row.names(sorted.table)
#'
#' # hierarchical clustering, and cutting into 5 clusters
#' hc <- hclust_haystack(dat.tsne, detection=dat.detection, genes=gene.subset)
#' hc.clusters <- cutree(hc,k = 5)
hclust_haystack_raw = function(x, y, detection, genes, method="ward.D"){

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(!is.numeric(y))
    stop("'y' must be a numeric vector")
  if(length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if(ncol(detection) != length(x))
    stop("The number of columns in 'detection' must be the same as the length of 'x'")
  if(!all(is.character(genes)))
    stop("Value of 'genes' should be characters")
  if(sum(is.element(genes, rownames(detection)))==0)
    stop("None of the values in 'genes' are present in row names of 'detection'")

  # if detection is a lgCMatrix, convert it to a lgRMatrix
  if(inherits(detection, "lgCMatrix")){
    message("### converting detection data from lgCMatrix to lgRMatrix")
    detection <- as(detection, "RsparseMatrix")
  }

  # get densities (not in high resolution)
  detection.rownames <- rownames(detection)
  row.index.subset <- which(is.element(detection.rownames, genes))

  densities <- get_density(x=x, y=y, detection=detection, rows.subset = row.index.subset)
  mat.dens <- apply(densities,1, function(x) as.vector(x))
  colnames(mat.dens) <- detection.rownames[row.index.subset]

  dist <- as.dist(1 - cor(mat.dens))
  hc <- hclust(dist, method=method)
  hc
}




### kmeans_haystack
#' Function for k-means clustering of genes according to their distribution on a 2D plot.
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param detection A logical matrix showing which genes (rows) are detected in which cells (columns)
#' @param genes A set of genes (of the 'detection' data) which will be clustered.
#' @param k The number of clusters to return.
#' @param ... Additional parameters which will be passed on to the kmeans function.
#'
#' @return An object of class kmeans, describing a clustering into 'k' clusters
#' @export
#'
#' @examples
#' # using the toy example of the singleCellHaystack package
#' # define a logical matrix with detection of each gene (rows) in each cell (columns)
#' dat.detection <- dat.expression > 1
#'
#' # running haystack in default mode
#' res <- haystack(dat.tsne, detection=dat.detection, method = "2D")
#'
#' # get biased genes, store in variable gene.subset
#' sorted.table <- show_result_haystack(res, p.value.threshold = 1e-5)
#' gene.subset <- row.names(sorted.table)
#'
#' # k-means clustering into 5 clusters
#' km <- kmeans_haystack(dat.tsne, detection=dat.detection, genes=gene.subset, k=5)
#' km.clusters <- km$cluster
kmeans_haystack_raw = function(x, y, detection, genes, k, ...){

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(!is.numeric(y))
    stop("'y' must be a numeric vector")
  if(length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if(ncol(detection) != length(x))
    stop("The number of columns in 'detection' must be the same as the length of 'x'")
  if(!all(is.character(genes)))
    stop("Value of 'genes' should be characters")
  if(sum(is.element(genes, rownames(detection)))==0)
    stop("None of the values in 'genes' are present in row names of 'detection'")
  if(missing(k) | !is.numeric(k) | k < 1)
    stop("Value of 'k' should be an integer larger than 1")

  # if detection is a lgCMatrix, convert it to a lgRMatrix
  if(inherits(detection, "lgCMatrix")){
    message("### converting detection data from lgCMatrix to lgRMatrix")
    detection <- as(detection, "RsparseMatrix")
  }

  # get densities (not in high resolution)
  detection.rownames <- rownames(detection)
  row.index.subset <- which(is.element(detection.rownames, genes))

  densities <- get_density(x=x, y=y, detection=detection, rows.subset = row.index.subset)
  mat.dens <- apply(densities,1, function(x) as.vector(x))
  colnames(mat.dens) <- detection.rownames[row.index.subset]

  km <- kmeans(x=t(mat.dens), centers=k, ...)
  km
}

