
########################################
########################################
### loading necessary packages

# for b splines
library(splines)


########################################
########################################
### kde2d_faster
#' Based on the MASS kde2d() function, but simplified beyond recognition; it's just tcrossprod() now.
#'
#' @param dens.x Contribution of all cells to densities of the x-axis grid points.
#' @param dens.y Contribution of all cells to densities of the y-axis grid points.
#'
kde2d_faster = function (dens.x, dens.y){
  tcrossprod(dens.x, dens.y)
}


########################################
########################################
#' Default function given by function bandwidth.nrd in MASS. No changes were made to this function.
#'
#' @param x A numeric vector
#'
#' @return A suitable bandwith.
default_bandwidth.nrd = function(x){
  r <- quantile(x, c(0.25, 0.75))
  h <- (r[2] - r[1])/1.34
  4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
}


########################################
########################################
#' Function that decides most of the parameters that will be during the "Haystack" analysis.
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param high.resolution Logical: should high resolution be used? Default is FALSE.
#'
#' @return A list containing various parameters to use in the analysis.
get_parameters_haystack = function(x,y,high.resolution=F){

  # I want to have a fixed number of bins between the 10th and 90th percentile in each dimension
  # this is to avoid bins getting squished because of a few outliers

  # the default number of bins between the 10% and 90% points in each dimension
  if(high.resolution==T){
    grid.points.10.90 <- 125
  } else {
    grid.points.10.90 <- 25
  }

  # get the 10% and 90% points in each dimension
  x10 <- as.numeric(quantile(x,0.10))
  x90 <- as.numeric(quantile(x,0.90))
  y10 <- as.numeric(quantile(y,0.10))
  y90 <- as.numeric(quantile(y,0.90))

  # get the minimum and maximum of each dimension
  x.min <- min(x)
  x.max <- max(x)
  y.min <- min(y)
  y.max <- max(y)

  # set the bin size and total number of bins in each dimension
  x.bin.size <- (x90-x10)/grid.points.10.90
  y.bin.size <- (y90-y10)/grid.points.10.90
  x.bins.10.to.min <- ceiling(abs(x.min-x10)/x.bin.size)
  x.bins.90.to.max <- ceiling(abs(x.max-x90)/x.bin.size)
  y.bins.10.to.min <- ceiling(abs(y.min-y10)/y.bin.size)
  y.bins.90.to.max <- ceiling(abs(y.max-y90)/y.bin.size)

  # get limits, grid points, bandwidths
  x.lim.min <- x10 - (x.bin.size*x.bins.10.to.min)
  x.lim.max <- x90 + (x.bin.size*x.bins.90.to.max)
  y.lim.min <- y10 - (y.bin.size*y.bins.10.to.min)
  y.lim.max <- y90 + (y.bin.size*y.bins.90.to.max)

  # the total number of grid points for both dimensions
  x.total.grid.points <- round((x.lim.max-x.lim.min)/x.bin.size)
  y.total.grid.points <- round((y.lim.max-y.lim.min)/y.bin.size)
  grid.points <- c(x.total.grid.points,y.total.grid.points)

  # the limits of the grid in both dimensions
  limits <- c(x.lim.min,x.lim.max,y.lim.min,y.lim.max)

  # the bandwidths used for the kernel function
  bandwidth.x <- default_bandwidth.nrd(x)
  bandwidth.y <- default_bandwidth.nrd(y)
  bandwidths <- c(bandwidth.x,bandwidth.y)
  bandwidths <- bandwidths/4

  # getting the distances (in units of bandwidths) between all cells and all grid points
  gx <- seq.int(x.lim.min, x.lim.max, length.out = x.total.grid.points)
  gy <- seq.int(y.lim.min, y.lim.max, length.out = y.total.grid.points)
  ax <- outer(gx, x, "-")/bandwidths[1L]
  ay <- outer(gy, y, "-")/bandwidths[2L]

  # densities based on a simplified Gaussian
  dens.x <- matrix(exp(-0.5 * ax * ax), , length(x))
  dens.y <- matrix(exp(-0.5 * ay * ay), , length(x)) # length(x) should be the same as length(y), of course

  # return results
  list(
    grid.points = grid.points,
    limits      = limits,
    bandwidths  = bandwidths,
    dens.x      = dens.x,
    dens.y      = dens.y
  )
}


########################################
########################################
#' Calculates the Kullback-Leibler divergence between distributions.
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param classes A logical vector. Values are T is the gene is expressed in a cell, F is not.
#' @param parameters Parameters of the analysis, as set by function 'get_parameters_haystack'
#' @param reference.prob A reference distribution to calculate the divergence against.
#' @param pseudo A pseudocount, used to avoid log(0) problems.
#'
#' @return A numerical value, the Kullback-Leibler divergence
get_D_KL = function(x, y, classes, parameters, reference.prob, pseudo){

  class.types = c(F,T)

  # the reference distribution Q of cells in the 2D plot
  Q <- reference.prob

  # calculating the Kullback-Leibler divergence of the distribution
  # of cells expressin and not expressing gene X vs reference distribution Q
  D_KLs <- c()
  for(c in 1:length(class.types)){
    cl <- class.types[c]
    cl.subset <- classes==cl
    if(sum(cl.subset)==0)
      next
    density <- kde2d_faster(dens.x=parameters$dens.x[,cl.subset],dens.y=parameters$dens.y[,cl.subset])

    # only decide pseudo if no value was given as input
    # (in total this takes some time)
    if(missing(pseudo))
      pseudo <- quantile(density[density>0],0.01)

    density <- density + pseudo
    P <- density / sum(density)

    D_KL <- sum(P * log(P/Q))
    D_KLs[c] <- D_KL
  }
  sum(D_KLs)
}


########################################
########################################
#' Estimates the significance of the observed Kullback-Leibler divergence by comparig to randomizations.
#'
#' @param T.counts The number of cells in which a gene is detected.
#' @param D_KL.observed A vector of observed Kullback-Leibler divergences.
#' @param D_KL.randomized A matrix of Kullback-Leibler divergences of randomized datasets.
#' @param output.dir Optional parameter. Default is NULL. If not NULL, some files will be written to this directory.
#'
#' @return A vector of log10 p values, corrected for multiple testing using the Bonferroni correction.
get_log_p_D_KL = function(T.counts, D_KL.observed, D_KL.randomized, output.dir = NULL){

  t.points <- as.numeric(row.names(D_KL.randomized))
  dat.mean.log2 <- apply(log2(D_KL.randomized),1,mean)
  dat.sd.log2 <- apply(log2(D_KL.randomized),1,sd)

  ##############################################
  ##############################################
  # fitting a spline for the mean D_KL

  # set df to 10, but lower if there are few points
  # however, df should be at least 5
  df.mean <- min(10, round(length(t.points)/20))
  df.mean <- max(df.mean,5)

  # set boundary knots to be slightly outside range of predictor values
  boundary.knots <- c(min(T.counts)-5, max(T.counts)+5)
  model.mean.log2 <- lm(dat.mean.log2 ~ bs(t.points, df=df.mean,Boundary.knots=boundary.knots))
  #summary(model.mean.log2)

  t.range <- range(t.points)
  t.seq <- seq(t.range[1],t.range[2],length.out=1000)
  fitted.mean.log2 <- predict(model.mean.log2, data.frame(t.points=t.seq), type="response")

  if(!is.null(output.dir)){
    outfile <- paste0(output.dir,"/fit.mean.log.D_KL_df",df.mean,".pdf")
    pdf(outfile)
    plot(t.points, dat.mean.log2)
    points(t.seq,fitted.mean.log2,col="red", type="l")
    dev.off()
  }


  ##############################################
  ##############################################
  # fitting a spline for the standard deviation of D_KL

  # set df to 10, but lower if there are few points
  # however, df should be at least 5
  df.sd <- min(10, round(length(t.points)/20))
  df.sd <- max(df.sd,5)

  model.sd.log2 <- lm(dat.sd.log2 ~ bs(t.points, df=df.sd,Boundary.knots=boundary.knots))
  #summary(model.sd.log2)

  fitted.sd.log2 <- predict(model.sd.log2, data.frame(t.points=t.seq), type="response")

  if(!is.null(output.dir)){
    outfile <- paste0(output.dir,"/fit.sd.log.D_KL_df",df.sd,".pdf")
    pdf(outfile)
    plot(t.points, dat.sd.log2)
    points(t.seq,fitted.sd.log2,col="red", type="l")
    dev.off()
  }

  ##############################################
  ##############################################
  # apply on actual observed values

  fitted.mean.log2 <- predict(model.mean.log2, data.frame(t.points=T.counts), type="response")
  fitted.sd.log2 <- predict(model.sd.log2, data.frame(t.points=T.counts), type="response")

  fitted.log.p.vals <- pnorm(log2(D_KL.observed), mean = fitted.mean.log2, sd = fitted.sd.log2, lower.tail = F, log.p = T)/log(10)

  # bonferroni correction for multiple testing
  fitted.log.p.vals <- fitted.log.p.vals + log10(length(fitted.log.p.vals))
  fitted.log.p.vals[fitted.log.p.vals>0] <- 0 # p values should be at most 1; so log10 should be <= 0

  fitted.log.p.vals
}


########################################
########################################
#' The main Haystack function.
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param logical A logical matrix showing which gens (rows) are detected in which cells (columns)
#' @param use.advanced.sampling If NULL naive sampling is used. If a vector is given (of length = no. of cells) sampling is done according to the values in the vector.
#' @param dir.randomization If NULL, no output is made about the random sampling step. If not NULL, files related to the randomizations are printed to this directory.
#'
#' @return An object of class "haystack"
#' @export
#'
#' @examples
#' warn("I will add this later")
haystack = function(x, y, logical, use.advanced.sampling=NULL, dir.randomization = NULL){

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(!is.numeric(y))
    stop("'y' must be a numeric vector")
  if(length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if(ncol(logical) != length(x))
    stop("The number of columns in 'logical' must be the same as the length of 'x'")
  if(!is.null(use.advanced.sampling)){
    if(!is.numeric(use.advanced.sampling))
      stop("'use.advanced.sampling' must either be NULL or a numeric vector")
    if(length(use.advanced.sampling) != length(x))
      stop("The length of 'use.advanced.sampling' must be the same as that of 'x'")
  }

  # warn about unusal input sizes
  if(length(x) < 50)
    warn("The number of cells seems very low (",length(x),"). Check your input.")
  if(nrow(logical) < 100)
    warn("The number of genes seems very low (",nrow(logical),"). Check your input.")


  # make dir if needed
  if(!is.null(dir.randomization)){
    if(!dir.exists(dir.randomization))
      dir.create(dir.randomization)
  }

  count.cells <- ncol(logical)
  count.genes <- nrow(logical)

  ##################################
  ### get reference probabilities "Q"
  ##################################

  # using all points, get number of grid points, limits, and bandwidths
  # get 2D densities using all points (using above grid points, limits, bandwidths)
  # add pseudocount to densities to avoid Inf problems
  # normalize to sum to 1
  message("### setting parameters...")
  parameters <- get_parameters_haystack(x,y)
  if(is.null(use.advanced.sampling)){
    density <- kde2d_faster(dens.x=parameters$dens.x,dens.y=parameters$dens.y)
    pseudo <- quantile(density[density>0],0.01)
  } else {
    density <- kde2d_faster(dens.x=t(t(parameters$dens.x)*use.advanced.sampling),
                            dens.y=t(t(parameters$dens.y)*use.advanced.sampling))
    density <- density / sum(density)
    pseudo <- quantile(density[density>0],0.01)
  }
  density <- density + pseudo
  # heatmap(density, Rowv=NA,Colv=NA,scale="none")
  Q <- density / sum(density)


  ##################################
  ### get probabilities "P" for each group ("F" and "T")
  ### this has to be one for every gene X
  ##################################

  # for class "F" points, and for "T" points, separately
  # get 2D densities (using above grid points, limits, bandwidths)
  # add pseudocount to densities to avoid Inf problems
  # normalize to sum to 1
  # get D_KL (or relative entropy) of this P vs reference Q
  message("### calculating Kulback-Leibler divergences...")
  D_KL.observed <- c()
  for(i in 1:count.genes){
    D_KL.observed[i] <- get_D_KL(x=x, y=y, classes=logical[i,], parameters=parameters, reference.prob=Q, pseudo=pseudo)
    if(i%%1000==0)
      message(paste0("### ... ",i," values out of ",count.genes," done"))
  }
  # return the sum of D_KL for "F" and "T"
  # store this value for each gene X


  ##################################
  ### use randomization to estimate expected values
  ### and p values for D_KL
  ##################################

  # Randomized D_KL values depend on the number of "F" and "T" cases
  # Thereofre, for all observed numbers of "T" cases (from 0 to 100%):
  # - do x randomizations and get teir D_KL
  # - get mean and SD per fraction,
  # - use those to estimate p values

  message("### starting randomizations...")
  T.counts <- apply(logical,1,sum)
  T.counts.unique <- sort(unique(T.counts))
  T.counts.unique.no <- length(T.counts.unique)
  p.vals <- rep(NA,nrow(logical))
  p.types <- rep(NA,nrow(logical))
  randomization.count <- 50

  # select T counts to include in randomizations
  T.counts.to.select <- 200
  if(T.counts.unique.no > T.counts.to.select){
    # random selection:
    # T.counts.selected <- sample(T.counts.unique, size=T.counts.to.select, replace = F)

    # more evenly spread selection:
    #tmp.T <- T.counts.unique[!T.counts.unique %in% c(0,count.cells)] # remove 0 or all, if present
    #tmp.T.no <- length(tmp.T)
    #T.counts.selected <- tmp.T[as.integer(seq(1,tmp.T.no,length.out = T.counts.to.select))]

    # include 10 lowest and highest values, otherwise evenly spread selection:
    tmp.T <- T.counts.unique[!T.counts.unique %in% c(0,count.cells)] # remove 0 or all, if present
    tmp.T.no <- length(tmp.T)
    T.counts.selected <- tmp.T[c(1:9,
                               as.integer(seq(10,tmp.T.no-9,length.out = T.counts.to.select-18)),
                               (tmp.T.no-8):tmp.T.no)
                               ]


  } else {
    # include all
    tmp.T <- T.counts.unique[!T.counts.unique %in% c(0,count.cells)] # remove 0 or all, if present
    tmp.T.no <- length(tmp.T)
    T.counts.selected <- tmp.T
    T.counts.to.select <- tmp.T.no
  }

  # to store randomization results
  all.D_KL.randomized <- matrix(NA,nrow=T.counts.to.select,ncol=randomization.count,
                                dimnames=list(T.counts.selected,1:randomization.count))

  # taking if - else statements outside of the for loop
  if(!is.null(use.advanced.sampling)){
    sampling.probs <- use.advanced.sampling

    for(i in 1:T.counts.to.select){
      if(i%%10==0)
        message(paste0("### ... ",i," values out of ",T.counts.to.select," done"))

      T.count <- T.counts.selected[i]

      D_KL.randomized <- rep(NA,randomization.count)
      for(r in 1:randomization.count){
        # using sampling according to the number of genes expressed in each cell
        # pick cells according to the number of genes they express
        samp <- sample(x=count.cells, prob=sampling.probs, size=T.count, replace = F)
        # turn into T or F
        classes <- is.element(1:count.cells,samp)
        D_KL.randomized[r] <- get_D_KL(x=x, y=y,
                                       classes=classes,
                                       parameters=parameters, reference.prob=Q, pseudo=pseudo)
      }
      all.D_KL.randomized[i,] <- D_KL.randomized
    }# end for all T counts to select
  } else {

    for(i in 1:T.counts.to.select){
      if(i%%10==0)
        message(paste0("### ... ",i," values out of ",T.counts.to.select," done"))

      T.count <- T.counts.selected[i]

      D_KL.randomized <- rep(NA,randomization.count)
      vector.to.randomize <- c(rep(T,T.count),rep(F,ncol(logical)-T.count))
      for(r in 1:randomization.count){
        # using default sampling
        D_KL.randomized[r] <- get_D_KL(x=x, y=y,
                                       classes=sample(x = vector.to.randomize),
                                       parameters=parameters, reference.prob=Q, pseudo=pseudo)
      }
      all.D_KL.randomized[i,] <- D_KL.randomized
    }# end for all T counts to select
  }# end if else

  message("### estimating p-values...")
  p.vals <- get_log_p_D_KL(T.counts = T.counts, D_KL.observed = D_KL.observed, D_KL.randomized = all.D_KL.randomized, output.dir = dir.randomization)

  if(!is.null(dir.randomization)){
    message("### writing randomized Kulback-Leibler divergences to file...")
    outputfile.randomized.D_KL <- paste0(dir.randomization,"/random.D_KL.csv")
    write.csv(file=outputfile.randomized.D_KL,all.D_KL.randomized)
  }

  message("### returning result...")
  # prepare the 'haystack' object to return
  res <- list(
    results = data.frame(
      D_KL = D_KL.observed,
      log.p.vals = p.vals,
      T.counts = T.counts,
      row.names = row.names(logical)
    )
  )
  class(res) <- "haystack"
  res
}


########################################
########################################
### get_density
# a function to get the density of points with value TRUE in the (x,y) plot
get_density = function(x, y, logical, rows.subset=1:nrow(logical), high.resolution=F){

  # set the parameters for getting the densities
  parameters <- get_parameters_haystack(x,y,high.resolution)

  densities <- array(data=NA, dim=c(length(rows.subset),parameters$grid.points))

  cl <- T # we are only looking at the T points here
  for(i in 1:length(rows.subset)){
    r <- rows.subset[i]
    x.subset <- logical[r,]==cl
    y.subset <- logical[r,]==cl
    density <- kde2d_faster(dens.x=parameters$dens.x[,x.subset],
                            dens.y=parameters$dens.y[,y.subset])

    densities[i,,] <- density

  }

  densities
}


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
