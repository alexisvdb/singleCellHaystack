

#' Visualizing the detection/expression of a gene in a 2D plot
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param gene name of a gene that is present in the input expression data, or a numerical index
#' @param expression a logical/numerical matrix showing detection/expression of genes (rows) in cells (columns)
#' @param detection an optional logical matrix showing detection of genes (rows) in cells (columns). If left as NULL, the density distribution of the gene is not plotted.
#' @param high.resolution logical (default: FALSE). If set to TRUE, the density plot will be of a higher resolution
#'
#' @return A plot
#' @export
#'
#' @examples
#' warning("I will add this later")
plot_gene_haystack = function(x, y, gene, expression, detection = NULL, high.resolution=F){

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
  if(!is.integer(gene) & !is.character(gene))
    stop("Value of 'gene' (\"",gene,"\") should be either a character, or an integer")
  if(is.integer(gene) & gene > nrow(expression))
    stop("Integer value of 'gene' (\"",gene,"\") is larger than the number of rows in 'expression' (",nrow(expression),")")
  if(is.character(gene) & !is.element(gene, rownames(expression)))
    stop("Value of 'gene' (string \"",gene,"\") is not present in row names of 'expression'")
  if(!is.logical(high.resolution))
    stop("Value of 'high.resolution' should be logical (TRUE or FALSE")

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

  # if expression is logical: treat as detection (TRUE or FALSE)
  # else, treat as expression levels
  if(is.numeric(expression)){
    Level <- expression[gene,]
    d <- d + geom_point(data=data.frame(x=x,y=y), aes(x, y, colour=Level)) + scale_color_gradient(low="grey", high="red")
  } else {
    Detection <- expression[gene,]
    d <- d + geom_point(data=data.frame(x=x,y=y), aes(x, y,colour=Detection))
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
#'
#' @return A plot
#' @export
#'
#' @examples
#' warning("I will add this later")
plot_gene_set_haystack = function(x, y, genes=NA, detection, high.resolution=T){

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(!is.numeric(y))
    stop("'y' must be a numeric vector")
  if(length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if(ncol(detection) != length(x))
    stop("The number of columns in 'detection' must be the same as the length of 'x'")
  if(!all(is.na(genes)) & !all(is.integer(genes)) & !all(is.character(genes)))
    stop("Entries of 'genes' should be either characters, or integers")
  if(all(is.integer(genes)) & any(genes > nrow(detection)))
    stop("Integer index of some genes is larger than the number of rows in 'expression' (",nrow(detection),")")
  if(!all(is.na(genes)) & all(is.character(genes)) & sum(is.element(genes, rownames(detection)))==0)
    stop("None of the entries in 'genes' are present in row names of 'detection'")
  if(!is.logical(high.resolution))
    stop("Value of 'high.resolution' should be logical (TRUE or FALSE")

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
  Detection <- mean.detection
  d <- d + geom_point(data=data.frame(x=x,y=y), aes(x, y,colour=Detection)) + scale_color_gradient(low="grey", high="red")

  d
}


