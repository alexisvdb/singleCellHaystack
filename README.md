# single-cell-haystack

## Introduction
'haystack' is a package for finding surprising needles (=genes) in haystacks (=2D representations of single cell transcriptome data). Single-cell RNA-seq (scRNA-seq) data is often represented in 2-dimentional plots (e.g. plots of two principal components, or t-SNE plots). 'haystack' can be used for finding genes that are expressed in subsets of cells that are non-randomly distributed in this 2D representation.

## Example usage

A small toy dataset is included in the package. The toy dataset includes:

- 'dat.expression': scRNA-seq expression of genes (rows) in cells (columns)

- 'dat.tsne':       a 2D representation of the cell in dat.expression


``` r
# Turn the expression data into detection (gene detected = T, not detected = F)
dat.detection <- dat.expression > 1

# run the main 'haystack' analysis
res <- haystack(x=dat.tsne$tSNE1, y=dat.tsne$tSNE2, detection=dat.detection)

# the returned results 'res' is of class 'haystack'
class(res)
## [1] "haystack"

# show top 10 "surprising" genes
show_result_haystack(res.haystack = res, n=10)

# alternatively: use a p-value threshold
#show_result_haystack(res.haystack = res, p.value.threshold = 1e-10)

# visualize one of the surprizing genes
plot_gene_haystack(x=dat.tsne$tSNE1, y=dat.tsne$tSNE2, expression=dat.expression, 
                      gene=gene, detection = dat.detection, high.resolution = T)

# get the top most significant genes, and cluster them by their distribution pattern in the 2D plot
sorted.table <- show_result_haystack(res.haystack = res, p.value.threshold = 1e-10)
gene.subset <- row.names(sorted.table)

# k-means clustering
km <- kmeans_haystack(x=dat.tsne$tSNE1, y=dat.tsne$tSNE2, detection=dat.detection, genes=gene.subset, k=5)
km.clusters <- km$cluster

# alternatively: hierarchical clustering
#hc <- hclust_haystack(x=dat.tsne$tSNE1, y=dat.tsne$tSNE2, detection=dat.detection, genes=gene.subset)
#hc.clusters <- cutree(hc,k = 5)

# visualize cluster distributions
plot_gene_set_haystack(x=dat.tsne$tSNE1, y=dat.tsne$tSNE2, detection=dat.detection, 
                          genes=names(km.clusters[km.clusters==1]))

```

