

#' Visualizing the detection/expression of a gene in a 2D plot
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param gene name of a gene that is present in the input expression data, or a numerical index
#' @param expression a logical/numerical matrix showing detection/expression of genes (rows) in cells (columns)
#' @param detection an optional logical matrix showing detection of genes (rows) in cells (columns). If left as NULL, the density distribution of the gene is not plotted.
#' @param high.resolution logical (default: FALSE). If set to TRUE, the density plot will be of a higher resolution
#' @param point.size numerical value to set size of points in plot. Default is 1.
#' @param order.by.signal If TRUE, cells with higher signal will be put on the foreground in the plot. Default is FALSE.
#'
#' @return A plot
#' @export
#'
#' @examples
#' # using the toy example of the singleCellHaystack package
#' # define a logical matrix with detection of each gene (rows) in each cell (columns)
#' dat.detection <- dat.expression > 1
#'
#' # running haystack in default mode
#' res <- haystack(dat.tsne, detection=dat.detection, method = "2D")
#' # list top 10 biased genes
#' show_result_haystack(res, n =10)
#'
#' # various was of plotting gene expression patterns
#' plot_gene_haystack(dat.tsne, expression=dat.expression, gene="gene_242",
#'  detection = dat.detection, high.resolution = TRUE)
#' plot_gene_haystack(dat.tsne, expression=dat.expression, gene="gene_242",
#'  detection = dat.detection, high.resolution = TRUE, point.size = .1)
#' plot_gene_haystack(dat.tsne, expression=dat.expression, gene="gene_242",
#'  high.resolution = TRUE)
#' plot_gene_haystack(dat.tsne, expression=dat.detection, gene="gene_242",
#'  detection = dat.detection, high.resolution = TRUE)
#' plot_gene_haystack(dat.tsne, expression=dat.detection, gene="gene_242",
#'  high.resolution = TRUE)
#'
#' # sort cells in the plot so cells with high signal come on top
#' plot_gene_haystack(dat.tsne, expression=dat.expression, gene="gene_242",
#'  high.resolution = TRUE, point.size=2, order.by.signal=TRUE)
plot_gene_haystack_raw = function(x, y, gene, expression, detection = NULL, high.resolution=FALSE, point.size=1, order.by.signal=FALSE){

  if(is.data.frame(expression)){
    warning("'expression' is a matrix, should be a matrix. Converting to matrix.")
    expression <- as.matrix(expression)
  }

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(!is.numeric(y))
    stop("'y' must be a numeric vector")
  if(length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if(ncol(expression) != length(x))
    stop("The number of columns in 'expression' must be the same as the length of 'x'")
  if(!is.null(detection)){
    if(!is.logical(detection))
      stop("'detection' must be either a logical matrix or NULL")
    if(any(dim(detection)!=dim(expression)))
      stop("'detection' must be of the same dimensions as 'expression', or should be NULL")
  }
  if(!is.numeric(gene) & !is.character(gene))
    stop("Value of 'gene' should be either a character, or an integer")
  if(is.numeric(gene) & gene > nrow(expression))
    stop("Integer value of 'gene' (\"",gene,"\") is larger than the number of rows in 'expression' (",nrow(expression),")")
  if(is.character(gene) & !is.element(gene, rownames(expression)))
    stop("Value of 'gene' (string \"",gene,"\") is not present in row names of 'expression'")
  if(!is.logical(high.resolution))
    stop("Value of 'high.resolution' should be logical (TRUE or FALSE")
  if(!is.numeric(point.size))
    stop("'point.size' must have a numeric value")
  if(!is.logical(order.by.signal))
    stop("Value of 'order.by.signal' should be logical (TRUE or FALSE")

  # set index of gene if it is a character
  # else, if it is an integer, use its value as an index
  if(is.character(gene)){
    gene.index <- which(rownames(expression)==gene)[1]
  } else {
    gene.index <- gene
  }

  # if detection is NULL: just show expression or presence of the gene in the 2D plot
  # else: get the density distribution of the cells expressing the gene and add it to the plot
  if(is.null(detection)){
    d <- ggplot() + theme_bw()
  } else {
    # get the density pattern around cells
    dens <- get_density(x=x, y=y, detection=detection, rows.subset=gene.index, high.resolution = high.resolution)
    dens.melted <- melt(dens[1,,])
    colnames(dens.melted) <- c("x", "y", "Density")
    d <- ggplot(dens.melted, aes(x, y))
    d <- d + geom_raster(aes_string(fill = "Density")) + scale_fill_gradient(low = "white", high = "steelblue")
    d <- d + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  }
  d <- d + theme(panel.border = element_rect(colour = "black", fill=NA))

  # prepare to plot "expression" levels
  dat.plot <- data.frame(
    x     = x,
    y     = y,
    value = unlist(expression[gene,])
  )

  # if needed order by signal
  if(order.by.signal){
    o <- order(dat.plot$value,decreasing=F)
    dat.plot <- dat.plot[o,]
  }

  # if expression is logical: treat as detection (TRUE or FALSE)
  # else, treat as expression levels
  if(is.numeric(dat.plot$value)){
    d <- d + geom_point(data=dat.plot, aes(x, y, colour=value), size=point.size) + scale_color_gradient(low="grey", high="red") + labs(color = "Level")
  } else if(is.logical(dat.plot$value)) {
    d <- d + geom_point(data=dat.plot, aes(x, y, colour=value), size=point.size) + labs(color = "Detection")
  }

  d
}







#' Visualizing the detection/expression of a set of genes in a 2D plot
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param genes Gene names that are present in the input expression data, or a numerical indeces. If NA, all genes will be used.
#' @param detection a logical matrix showing detection of genes (rows) in cells (columns)
#' @param high.resolution logical (default: TRUE). If set to FALSE, the density plot will be of a lower resolution
#' @param point.size numerical value to set size of points in plot. Default is 1.
#' @param order.by.signal If TRUE, cells with higher signal will be put on the foreground in the plot. Default is FALSE.
#'
#' @return A plot
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
#' hc <- hclust_haystack(dat.tsne, detection=dat.detection,
#'  genes=gene.subset)
#' hc.clusters <- cutree(hc,k = 5)
#'
#' # visualization of average pattern of cluster 1
#' plot_gene_set_haystack(dat.tsne, detection=dat.detection,
#'  genes=names(hc.clusters[hc.clusters==1]))
#'
#' # tweak size of points in plot sing 'point.size'
#' plot_gene_set_haystack(dat.tsne, detection=dat.detection,
#'  genes=names(hc.clusters[hc.clusters==1]), point.size=.1)
#'
#' # sort cells in the plot so cells with high average signal come on top
#' plot_gene_set_haystack(dat.tsne, detection=dat.detection,
#'  genes=names(hc.clusters[hc.clusters==1]), point.size=2,
#'  order.by.signal=TRUE)
plot_gene_set_haystack_raw = function(x, y, genes=NA, detection, high.resolution=TRUE, point.size=1, order.by.signal=FALSE){

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(!is.numeric(y))
    stop("'y' must be a numeric vector")
  if(length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if(ncol(detection) != length(x))
    stop("The number of columns in 'detection' must be the same as the length of 'x'")
  if(!all(is.na(genes)) & !all(is.numeric(genes)) & !all(is.character(genes)))
    stop("Entries of 'genes' should be either characters, or integers")
  if(all(is.numeric(genes)) & any(genes > nrow(detection)))
    stop("Integer index of some genes is larger than the number of rows in 'expression' (",nrow(detection),")")
  if(!all(is.na(genes)) & all(is.character(genes)) & sum(is.element(genes, rownames(detection)))==0)
    stop("None of the entries in 'genes' are present in row names of 'detection'")
  if(!is.logical(high.resolution))
    stop("Value of 'high.resolution' should be logical (TRUE or FALSE")
  if(!is.numeric(point.size))
    stop("'point.size' must have a numeric value")
  if(!is.logical(order.by.signal))
    stop("Value of 'order.by.signal' should be logical (TRUE or FALSE")

  # set index of gene if it is a character
  # else, if it is an integer, use its value as an index
  if(all(is.na(genes))){
    gene.indices <- 1:nrow(detection)
  }
  else if(all(is.character(genes))){
    gene.indices <- which(is.element(rownames(detection),genes))
  } else {
    gene.indices <- unique(genes)
  }


  ### get mean density

  # get densities of all genes
  dens <- get_density(x=x, y=y, detection=detection, rows.subset = gene.indices, high.resolution = high.resolution)

  # get average density
  mean.density <- apply(dens,c(2,3),mean)


  ### get mean detection level
  if(length(gene.indices)==1){
    mean.detection <- apply(t(detection[gene.indices,]),2,mean)
  } else {
    mean.detection <- apply(detection[gene.indices,],2,mean)
  }


  ### make plot
  dens.melted <- melt(mean.density)
  colnames(dens.melted) <- c("x", "y", "Density")
  d <- ggplot(dens.melted, aes(x, y)) +
    theme(panel.border = element_rect(colour = "black", fill=NA))
  d <- d + geom_raster(aes_string(fill = "Density")) + scale_fill_gradient(low = "white", high = "steelblue")
  d <- d + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

  # prepare to plot "expression" levels
  dat.plot <- data.frame(
    x     = x,
    y     = y,
    value = mean.detection
  )

  # if needed order by signal
  if(order.by.signal){
    o <- order(dat.plot$value,decreasing=F)
    dat.plot <- dat.plot[o,]
  }

  d <- d + geom_point(data=dat.plot, aes(x, y, colour=value), size=point.size) + scale_color_gradient(low="grey", high="red") + labs(color = "Detection")

  d
}


