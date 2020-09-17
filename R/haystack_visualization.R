

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
#' # various was of plotting gene expression patterns
#' plot_gene_haystack(dat.tsne, expression=dat.expression, gene="gene_242",
#'  detection = dat.detection, high.resolution = TRUE)
#' plot_gene_haystack(dat.tsne, expression=dat.expression, gene="gene_242",
#'  detection = dat.detection, high.resolution = TRUE, point.size = .1)
plot_gene_haystack_raw = function(x, y, gene, expression, detection = NULL, high.resolution = FALSE, point.size = 1, order.by.signal = FALSE){

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
    if(!is.matrix(detection) && ! inherits(detection, "lgCMatrix") && ! inherits(detection, "lgRMatrix"))
      stop("'detection' must be a matrix, lgCMatrix, or lgRMatrix")
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

  # if detection is a lgCMatrix, convert it to a lgRMatrix
  if(inherits(detection, "lgCMatrix")){
    message("### converting detection data from lgCMatrix to lgRMatrix")
    detection <- as(detection, "RsparseMatrix")
  }
  # if expression is a dgCMatrix, convert it to a dgRMatrix
  if(inherits(expression, "dgCMatrix") ){
    message("### converting expression data from dgCMatrix to dgRMatrix")
    expression <- as(expression, "RsparseMatrix")
  }
  # expression could also be logical
  # if expression is a lgCMatrix, convert it to a lgRMatrix
  if(inherits(expression, "lgCMatrix") ){
    message("### converting expression data from lgCMatrix to lgRMatrix")
    expression <- as(expression, "RsparseMatrix")
  }
  # expression could also be a dgeMatrix
  # if expression is a dgeMatrix, convert it to a dgRMatrix
  if(inherits(expression, "dgeMatrix") ){
    message("### converting expression data from dgeMatrix to lgRMatrix")
    expression <- as(expression, "RsparseMatrix")
  }

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
    d <- ggplot(dens.melted, aes_string("x", "y"))
    d <- d + geom_raster(aes_string(fill = "Density")) + scale_fill_gradient(low = "white", high = "steelblue")
    d <- d + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  }
  d <- d + theme(panel.border = element_rect(colour = "black", fill=NA))

  # prepare to plot "expression" levels
  dat.plot <- data.frame(
    x     = x,
    y     = y,
    value = if(inherits(expression, "dgRMatrix")){ # for numeric sparse case
              extract_row_dgRMatrix(expression, gene.index)
            } else if(inherits(expression, "lgRMatrix")){ # for logical sparse case
              extract_row_lgRMatrix(expression, gene.index)
            } else {
              unlist(expression[gene.index,])
            }
  )

  # if needed order by signal
  if(order.by.signal){
    o <- order(dat.plot$value,decreasing = FALSE)
    dat.plot <- dat.plot[o,]
  }

  # if expression is logical: treat as detection (TRUE or FALSE)
  # else, treat as expression levels
  if(is.numeric(dat.plot$value)){
    # d <- d + geom_point(data=dat.plot, aes_string("x", "y", colour="value"), size=point.size) + scale_color_gradient(low="grey", high="red") + labs(color = "Level")
    d <- d + geom_point(data=dat.plot, aes_string("x", "y", colour="value"), size=point.size) + scale_color_gradientn(colours=c("#BBBBBB","#F0E442","#E69F00","#FF0000")) + labs(color = "Level")
  } else if(is.logical(dat.plot$value)) {

    d <- d + geom_point(data=dat.plot, aes_string("x", "y", colour="value"), size=point.size) + labs(color = "Detection")
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
#' # define a set of genes that we want to visualize
#' # this might be a set of differnentially expressed genes
#' # predicted by haystack and clustered together by hclust_haystack
#' gene_set <- c("gene_9", "gene_59", "gene_112", "gene_137", "gene_155",
#'   "gene_216", "gene_234", "gene_275", "gene_291", "gene_317",
#'   "gene_339", "gene_340", "gene_351", "gene_400", "gene_424", "gene_479")
#'
#' # visualize the expression pattern of the set of genes
#' plot_gene_set_haystack(dat.tsne, detection=dat.detection, genes=gene_set)
plot_gene_set_haystack_raw = function(x, y, genes=NA, detection, high.resolution = TRUE, point.size=1, order.by.signal = FALSE){

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

  # if detection is a lgCMatrix, convert it to a lgRMatrix
  if(inherits(detection, "lgCMatrix")){
    message("### converting detection data from lgCMatrix to lgRMatrix")
    detection <- as(detection, "RsparseMatrix")
  }

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
    mean.detection <- Matrix::colMeans(detection[gene.indices,])
  }


  ### make plot
  dens.melted <- melt(mean.density)
  colnames(dens.melted) <- c("x", "y", "Density")
  d <- ggplot(dens.melted, aes_string("x", "y")) +
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
    o <- order(dat.plot$value,decreasing = FALSE)
    dat.plot <- dat.plot[o,]
  }

  #d <- d + geom_point(data=dat.plot, aes_string("x", "y", colour="value"), size=point.size) + scale_color_gradient(low="grey", high="red") + labs(color = "Detection")
  #d <- d + geom_point(data=dat.plot, aes_string("x", "y", colour="value"), size=point.size) + scale_color_gradientn(colours=c("grey","red","dark red")) + labs(color = "Detection")
  d <- d + geom_point(data=dat.plot, aes_string("x", "y", colour="value"), size=point.size) + scale_color_gradientn(colours=c("#BBBBBB","#F0E442","#E69F00","#FF0000")) + labs(color = "Detection")

  d
}


