

#' Function for hierarchical clustering of genes according to their distribution in a higher-dimensional space.
#'
#' @param x Coordinates of cells in a 2D or higher-dimensional space. Rows represent cells, columns the dimensions of the space.
#' @param detection A logical matrix showing which genes (rows) are detected in which cells (columns)
#' @param genes A set of genes (of the 'detection' data) which will be clustered.
#' @param method The method to use for hierarchical clustering. See '?hclust' for more information. Default: "ward.D".
#' @param grid.coordinates Coordinates of grid points in the same space as 'x', to be used to estimate densities for clustering.
#' @param scale whether to scale data.
#'
#' @return An object of class hclust, describing a hierarchical clustering tree.
#' @export
#'
#' @examples
#' # to be added
hclust_haystack_highD = function(x, detection, genes, method="ward.D", grid.coordinates = NULL, scale = TRUE){

  # if data.frame, convert to matrix
  if(is.data.frame(x))
    x <- as.matrix(x)

  # check input
  if(!is.numeric(x) | !is.matrix(x))
    stop("'x' must be a numeric matrix")
  if(ncol(x) < 2)
    stop("'x' must have at least 2 columns")
  if(ncol(detection) != nrow(x))
    stop("The number of columns in 'detection' must be the same as the rows in 'x'")
  if(!all(is.character(genes)))
    stop("Value of 'genes' should be characters")
  if(sum(is.element(genes, rownames(detection)))==0)
    stop("None of the values in 'genes' are present in row names of 'detection'")
  if(is.null(grid.coordinates))
    stop("Please provide 'grid.coordinates'. For example, the result of haystack_highD() includes such coordinates. ")
  if(ncol(grid.coordinates)!=ncol(x))
    stop("The number of columns in 'x' and 'grid.coordinates' don't match.")

  if(!is.logical(scale) | length(scale) > 1)
    stop("The value of 'scale' must be either TRUE or FALSE")

  # scale data if needed
  if(scale) {
    x <- scale(x)
    # save the mean and stdev of the scaling
    x.scale.center <- attr(x = x, which = "scaled:center")
    x.scale.scale <- attr(x = x, which = "scaled:scale")
    grid.coordinates <- (grid.coordinates - rep(x.scale.center,each=nrow(grid.coordinates))) / rep(x.scale.scale,each=nrow(grid.coordinates))
  }

  # if detection is a lgCMatrix, convert it to a lgRMatrix
  if(inherits(detection, "lgCMatrix")){
    message("### converting detection data from lgCMatrix to lgRMatrix")
    detection <- as(detection, "RsparseMatrix")
  }

  # get densities
  detection.rownames <- rownames(detection)
  row.index.subset <- which(is.element(detection.rownames, genes))

  dist.to.grid <- get_dist_two_sets(x,grid.coordinates)

  # process the distances to a suitable density contribution
  # first, set bandwidth
  # bandwidth <- sqrt(sum((apply(x, 2, default_bandwidth.nrd)) ^ 2))
  bandwidth <- median(apply(dist.to.grid,1,min))
  dist.to.grid.norm <- dist.to.grid / bandwidth
  density.contributions <-
    exp(-dist.to.grid.norm * dist.to.grid.norm / 2)

  densities <- matrix(NA, nrow=length(row.index.subset), ncol=ncol(density.contributions))
  row.names(densities) <- detection.rownames[row.index.subset]

  message("### collecting density data...")
  pb <- txtProgressBar(min = 0, max = length(row.index.subset), style = 3, file = stderr()) # progress bar
  if(is.matrix(detection)){
    for(g in 1:length(row.index.subset)){
      gene_index <- row.index.subset[g]
      densities[g,] <- apply(density.contributions[detection[gene_index,],],2,sum)
      setTxtProgressBar(pb, g) # progress bar
    }
  } else if( inherits(detection, "lgRMatrix") ){
    for(g in 1:length(row.index.subset)){
      gene_index <- row.index.subset[g]
      densities[g,] <- apply(density.contributions[extract_row_lgRMatrix(detection,gene_index),],2,sum)
      setTxtProgressBar(pb, g) # progress bar
    }
  } else {
    stop("'detection' must be a matrix or lgRMatrix")
  }
  close(pb) # progress bar

  #heatmap(dist.to.grid.norm, Rowv=NA, Colv=NA, scale="none")
  #heatmap(densities, Rowv=NA, Colv=NA, scale="none")

  # rescale to sum to 1. This is to avoid R thinking sd=0 in the case where an entire row has very low values
  densities <- densities / apply(densities,1,sum)

  dist <- as.dist(1 - cor(t(densities),method = "spearman")) # dist(densities)
  hc <- hclust(dist, method=method)
  hc
}



#' Function for k-means clustering of genes according to their distribution in a higher-dimensional space.
#'
#' @param x Coordinates of cells in a 2D or higher-dimensional space. Rows represent cells, columns the dimensions of the space.
#' @param detection A logical matrix showing which genes (rows) are detected in which cells (columns)
#' @param genes A set of genes (of the 'detection' data) which will be clustered.
#' @param grid.coordinates Coordinates of grid points in the same space as 'x', to be used to estimate densities for clustering.
#' @param k The number of clusters to return.
#' @param scale whether to scale data.
#' @param ... Additional parameters which will be passed on to the kmeans function.
#'
#' @return An object of class kmeans, describing a clustering into 'k' clusters
#' @export
#'
#' @examples
#' # to be added
kmeans_haystack_highD = function(x, detection, genes, grid.coordinates = NULL, k, scale = TRUE, ...){

  # if data.frame, convert to matrix
  if(is.data.frame(x))
    x <- as.matrix(x)

  # check input
  if(!is.numeric(x) | !is.matrix(x))
    stop("'x' must be a numeric matrix")
  if(ncol(x) < 2)
    stop("'x' must have at least 2 columns")
  if(ncol(detection) != nrow(x))
    stop("The number of columns in 'detection' must be the same as the rows in 'x'")
  if(!all(is.character(genes)))
    stop("Value of 'genes' should be characters")
  if(sum(is.element(genes, rownames(detection)))==0)
    stop("None of the values in 'genes' are present in row names of 'detection'")
  if(missing(k) | !is.numeric(k) | k < 1)
    stop("Value of 'k' should be an integer larger than 1")
  if(is.null(grid.coordinates))
    stop("Please provide 'grid.coordinates'. For example, the result of haystack_highD() includes such coordinates. ")
  if(ncol(grid.coordinates)!=ncol(x))
    stop("The number of columns in 'x' and 'grid.coordinates' don't match.")

  if(!is.logical(scale) | length(scale) > 1)
    stop("The value of 'scale' must be either TRUE or FALSE")

  # scale data if needed
  if(scale){
    x <- scale(x)
    # save the mean and stdev of the scaling
    x.scale.center <- attr(x = x, which = "scaled:center")
    x.scale.scale <- attr(x = x, which = "scaled:scale")
    grid.coordinates <- (grid.coordinates - rep(x.scale.center,each=nrow(grid.coordinates))) / rep(x.scale.scale,each=nrow(grid.coordinates))
  }

  # if detection is a lgCMatrix, convert it to a lgRMatrix
  if(inherits(detection, "lgCMatrix")){
    message("### converting detection data from lgCMatrix to lgRMatrix")
    detection <- as(detection, "RsparseMatrix")
  }

  # get densities
  detection.rownames <- rownames(detection)
  row.index.subset <- which(is.element(detection.rownames, genes))

  dist.to.grid <- get_dist_two_sets(x,grid.coordinates)

  # process the distances to a suitable density contribution
  # first, set bandwidth
  # bandwidth <- sqrt(sum((apply(x, 2, default_bandwidth.nrd)) ^ 2))
  bandwidth <- median(apply(dist.to.grid,1,min))
  dist.to.grid.norm <- dist.to.grid / bandwidth
  density.contributions <-
    exp(-dist.to.grid.norm * dist.to.grid.norm / 2)

  densities <- matrix(NA, nrow=length(row.index.subset), ncol=ncol(density.contributions))
  row.names(densities) <- detection.rownames[row.index.subset]

  message("### collecting density data...")
  pb <- txtProgressBar(min = 0, max = length(row.index.subset), style = 3, file = stderr()) # progress bar
  if(is.matrix(detection)){
    for(g in 1:length(row.index.subset)){
      gene_index <- row.index.subset[g]
      densities[g,] <- apply(density.contributions[detection[gene_index,],],2,sum)
      setTxtProgressBar(pb, g) # progress bar
    }
  } else if( inherits(detection, "lgRMatrix") ){
    for(g in 1:length(row.index.subset)){
      gene_index <- row.index.subset[g]
      densities[g,] <- apply(density.contributions[extract_row_lgRMatrix(detection,gene_index),],2,sum)
      setTxtProgressBar(pb, g) # progress bar
    }
  } else {
    stop("'detection' must be a matrix or lgRMatrix")
  }
  close(pb) # progress bar

  km <- kmeans(x=densities, centers=k, ...)
  km
}



