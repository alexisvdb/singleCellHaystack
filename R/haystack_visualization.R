
########################################
########################################
### loading necessary packages

# for plots
library(ggplot2)
library(reshape2)


########################################
########################################
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
#' warn("I will add this later")
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
    d <- ggplot()
  } else {
    # get the density pattern around cells
    dens <- get_density(x=dat.tsne$tSNE1, y=dat.tsne$tSNE2, logical=detection, rows.subset=gene.index, high.resolution = high.resolution)
    dens.melted <- melt(dens[1,,])
    d <- ggplot(dens.melted, aes(Var1, Var2))
    d <- d + geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "steelblue")
  }

  # if expression is logical: treat as detection (TRUE or FALSE)
  # else, treat as expression levels
  if(is.numeric(expression)){
    Level <- expression[gene,]
    d <- d + geom_point(data=dat.tsne, aes(tSNE1, tSNE2,colour=Level)) + scale_color_gradient(low="grey", high="red")
  } else {
    Detection <- expression[gene,]
    d <- d + geom_point(data=dat.tsne, aes(tSNE1, tSNE2,colour=Detection))
  }

  d
}





