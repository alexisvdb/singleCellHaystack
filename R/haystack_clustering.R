########################################
########################################
### hclust_haystack
#' Function for hierarchical clustering of genes according to their distribution on a 2D plot.
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param detection A logical matrix showing which gens (rows) are detected in which cells (columns)
#' @param genes A set of genes (of the 'detection' data) which will be clustered.
#' @param method The method to use for hierarchical clustering. See '?hclust' for more information. Default: "ward.D".
#'
#' @return An object of class hclust, describing a hierarchical clustering tree.
#' @export
#'
#' @examples
#' warning("I will add this later")
hclust_haystack= function(x, y, detection, genes, method="ward.D"){

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


  # get densities (not in high relosultion)
  detection.rownames <- rownames(detection)
  row.index.subset <- which(is.element(detection.rownames, genes))

  densities <- get_density(x=x, y=y, detection=detection, rows.subset = row.index.subset)
  mat.dens <- apply(densities,1, function(x) as.vector(x))
  colnames(mat.dens) <- detection.rownames[row.index.subset]

  dist <- as.dist(1 - cor(mat.dens))
  hc <- hclust(dist, method=method)
  hc
}


########################################
########################################
### kmeans_haystack
#' Function for k-means clustering of genes according to their distribution on a 2D plot.
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param detection A logical matrix showing which gens (rows) are detected in which cells (columns)
#' @param genes A set of genes (of the 'detection' data) which will be clustered.
#' @param k The number of clusters to return.
#' @param ... Additional parameters which will be passed on to the kmeans function.
#'
#' @return An object of class kmeans, describing a clustering into 'k' clusters
#' @export
#'
#' @examples
#' warning("I will add this later")
kmeans_haystack= function(x, y, detection, genes, k, ...){

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


  # get densities (not in high relosultion)
  detection.rownames <- rownames(detection)
  row.index.subset <- which(is.element(detection.rownames, genes))

  densities <- get_density(x=x, y=y, detection=detection, rows.subset = row.index.subset)
  mat.dens <- apply(densities,1, function(x) as.vector(x))
  colnames(mat.dens) <- detection.rownames[row.index.subset]

  km <- kmeans(x=t(mat.dens), centers=k, ...)
  km
}

