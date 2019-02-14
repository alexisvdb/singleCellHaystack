
########################################
########################################
### get_hierarchical_clustering
# clusters genes according to their density profile in the (x,y) plot
get_hierarchical_clustering = function(x, y, logical, rows.subset=1:nrow(logical)){
  densities <- get_density(x=x, y=y, logical=classes, rows.subset = rows.subset)
  mat.dens <- apply(densities,1, function(x) as.vector(x))
  colnames(mat.dens) <- row.names(logical)[rows.subset]
  dist <- as.dist(1 - cor(mat.dens))
  hc <- hclust(dist, method="ward.D")

  hc
}


########################################
########################################
### get_high_resolution_density_of_clusters
# Given clusters of genes, this function return the averaged density profile
# in the (x,y) plot for each cluster.
get_high_resolution_density_of_clusters = function(x, y, logical, clusters){

  unique.clusters <- sort(unique(clusters))
  mean.densities <- list()

  # run through the genes in each cluster and average their densities
  for(c in 1:length(unique.clusters)){
    cl <- unique.clusters[c]
    genes <- names(clusters[clusters==cl])
    gene.indices <- which(is.element(row.names(logical),genes))
    d <- get_density(x=x, y=y, logical=logical, rows.subset = gene.indices, high.resolution = T)
    mean.density <- apply(d,c(2,3),mean)
    mean.densities[[cl]] <- mean.density
    #heatmap.2(t(mean.density)[ncol(mean.density):1,], Rowv = NA, Colv=NA, dendrogram = "none", scale="none", trace="none")
  }

  mean.densities
}
