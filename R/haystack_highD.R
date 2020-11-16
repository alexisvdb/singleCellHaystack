#' Calculate the Euclidean distance between x and y.
#'
#' @param x A numerical vector.
#' @param y A numerical vector.
#'
#' @return A numerical value, the Euclidean distance.
get_euclidean_distance = function(x,y){
  sqrt(sum((x-y)^2))
}

#' Calculate the pairwise Euclidean distances between the rows of 2 matrices.
#'
#' @param set1 A numerical matrix.
#' @param set2 A numerical matrix.
#'
#' @return A matrix of pairwise distances between the rows of 2 matrices.
get_dist_two_sets = function(set1,set2){
  pairwise.dist <- matrix(NA,nrow=nrow(set1), ncol=nrow(set2))
  for(s2 in 1:nrow(set2)){
    point2 <- set2[s2,]
    pairwise.dist[,s2] <- apply(set1,1,function(x) get_euclidean_distance(x,point2))
  }# end for each point in set2
  pairwise.dist
}


#' Calculates the Kullback-Leibler divergence between distributions for the high-dimensional version of haystack().
#'
#' @param classes A logical vector. Values are T is the gene is expressed in a cell, F is not.
#' @param density.contributions A matrix of density contributions of each cell (rows) to each center point (columns).
#' @param reference.prob A reference distribution to calculate the divergence against.
#' @param pseudo A pseudocount, used to avoid log(0) problems.
#'
#' @return A numerical value, the Kullback-Leibler divergence
get_D_KL_highD = function(classes, density.contributions, reference.prob, pseudo = 0){

  class.types = c(FALSE, TRUE)

  # the reference distribution Q of cells in the 2D plot
  Q <- reference.prob

  # calculating the Kullback-Leibler divergence of the distribution
  # of cells expressing and not expressing gene X vs reference distribution Q
  D_KLs <- c()
  for(c in 1:length(class.types)){
    cl <- class.types[c]
    cl.subset <- classes==cl
    if(sum(cl.subset)==0){ # this means a gene is expressed or not expressed in all cells; return 0
      #D_KLs[c] <- 0
      return(0)
    } else {
      P <- apply(density.contributions[cl.subset, , drop = FALSE], 2, sum)
      P <- P + pseudo
      P <- P / sum(P)
      D_KL <- sum(P * log(P/Q))
      D_KLs[c] <- D_KL
    }
  }
  sum(D_KLs)
}



#' A function to decide grid points in a higher-dimensional space
#'
#' @param input A numerical matrix with higher-dimensional coordinates (columns) of points (rows)
#' @param method The method to decide grid points. Should be "centroid" (default) or "seeding".
#' @param grid.points The number of grid points to return. Default is 100.
#'
#' @return Coordinates of grid points in the higher-dimensonal space.
get_grid_points = function(input, method="centroid", grid.points = 100){

  if(nrow(input) < grid.points){
    warning("Fewer input points than grid.points. Adjusting value of grid.points to ",nrow(input))
    grid.points <- nrow(input)
  }

  if(method=="centroid"){
    # perform k-means clustering and get the centroids (centers) of each cluser
    # suppressing "did not converge" warnings, because convergence is not really important here
    suppressWarnings(
      res.kmeans <- kmeans(input, centers=grid.points, iter.max = 10, nstart = 10)
    )
    grid.coord <- res.kmeans$centers

  } else if(method=="seeding"){
    # do seeing as used in the seeding step of the kmeans++ algorithm
    # one by one pick points which are distal to the points picked so far

    input.dim <- ncol(input)
    input.points <- nrow(input)

    # pick a first grid point at random
    index <- sample(input.points,1)
    grid.coord <- matrix(NA, nrow=grid.points, ncol=input.dim)
    grid.coord[1,] <- input[index,]
    min.dist.to.grid <- get_dist_two_sets(set1 = input,set2 = grid.coord[1,,drop=F])

    # matrix with distances between inputs and grid points
    dist.to.grid <- matrix(NA, nrow=input.points, ncol = grid.points)
    dist.to.grid[,1] <- min.dist.to.grid

    for(i in 2:grid.points){

      # get probabilities based on min. distance to grid points
      probs <- min.dist.to.grid^2

      # pick new grid point
      index <- sample(input.points,1,prob = probs)
      grid.coord[i,] <- matrix(input[index,], nrow=1)

      # get distances to this new grid point and update min. distances
      tmp.dist.to.grid <- get_dist_two_sets(set1 = input,set2 = grid.coord[i,,drop=F])
      dist.to.grid[,i] <- tmp.dist.to.grid
      min.dist.to.grid <- apply(cbind(min.dist.to.grid,tmp.dist.to.grid),1,min)

    }

    # if necessary, average grid points out to nearest few points
    # here we always do this
    average.grid.points = TRUE
    if(average.grid.points){
      weights <- c(2,1,1) # more weight to closest point
      for(i in 1:grid.points){
        o <- order(dist.to.grid[,i])
        grid.coord[i,] <- apply(input[o[1:3],],2, weighted.mean, weights)
      }
    }

  }
  grid.coord

}




#' The main Haystack function, for higher-dimensional spaces.
#'
#' @param x Coordinates of cells in a 2D or higher-dimensional space. Rows represent cells, columns the dimensions of the space.
#' @param detection A logical matrix showing which genes (rows) are detected in which cells (columns)
#' @param use.advanced.sampling If NULL naive sampling is used. If a vector is given (of length = no. of cells) sampling is done according to the values in the vector.
#' @param dir.randomization If NULL, no output is made about the random sampling step. If not NULL, files related to the randomizations are printed to this directory.
#' @param scale Logical (default=TRUE) indicating whether input coordinates in x should be scaled to mean 0 and standard deviation 1.
#' @param grid.points An integer specifying the number of centers (gridpoints) to be used for estimating the density distributions of cells. Default is set to 100.
#' @param grid.method The method to decide grid points for estimating the density in the high-dimensional space. Should be "centroid" (default) or "seeding".
#'
#' @return An object of class "haystack", including the results of the analysis, and the coordinates of the grid points used to estimate densities.
#' @export
#'
#' @examples
#' # I need to add some examples.
#' # A toy example will be added too.
haystack_highD = function(x, detection, grid.points = 100, use.advanced.sampling=NULL, dir.randomization = NULL, scale=TRUE, grid.method="centroid"){
  message("### calling haystack_highD()...")

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric matrix")
  if(!is.matrix(x))
    stop("'x' must be a numeric matrix")
  if(ncol(x) < 2)
    stop("'x' must have at least 2 columns")
  if(!is.matrix(detection) && ! inherits(detection, "lgCMatrix") && ! inherits(detection, "lgRMatrix"))
    stop("'detection' must be a matrix, lgCMatrix, or lgRMatrix")
  if(ncol(detection) != nrow(x))
    stop("The number of columns in 'detection' must be the same as the rows in 'x'")
  if(!is.numeric(grid.points))
    stop("The value of 'grid.points' must be a numeric")
  if(grid.points >= ncol(detection))
    stop("The number of grid points appears to be very high (higher than the number of cells). You can set the number of grid points using the 'grid.points' parameter.")
  if(!is.null(use.advanced.sampling)){
    if(!is.numeric(use.advanced.sampling))
      stop("'use.advanced.sampling' must either be NULL or a numeric vector")
    if(length(use.advanced.sampling) != nrow(x))
      stop("The length of 'use.advanced.sampling' must be the same as the number of rows in 'x'")
  }
  if(!is.logical(scale) | length(scale) > 1)
    stop("The value of 'scale' must be either TRUE or FALSE")
  if(grid.method!="centroid" & grid.method!="seeding")
    stop("The value of 'grid.method' must be either 'centroid' or 'seeding'")

  # if detection is a lgCMatrix, convert it to a lgRMatrix
  if(inherits(detection, "lgCMatrix")){
    message("### converting detection data from lgCMatrix to lgRMatrix")
    detection <- as(detection, "RsparseMatrix")
  }

  count.cells <- ncol(detection)
  count.genes <- nrow(detection)

  # warn about unusual input sizes
  if(nrow(x) < 50)
    warning("The number of cells seems very low (",nrow(x),"). Check your input.")
  if(nrow(detection) < 100)
    warning("The number of genes seems very low (",nrow(detection),"). Check your input.")

  # warn about extreme values for 'grid.points'
  if(grid.points < 10)
    warning("The value of 'grid.points' appears to be very low (<10). You can set the number of grid points using the 'grid.points' parameter.")
  if(grid.points > count.cells/10)
    warning("The value of 'grid.points' appears to be very high (> No. of cells / 10). You can set the number of grid points using the 'grid.points' parameter.")

  # advanced sampling is slow on large datasets. Recommend using wrswoR
  if(!is.null(use.advanced.sampling)){
    if(requireNamespace("wrswoR", quietly = TRUE)){
      use.wrswoR <- TRUE
      message("### Using package wrswoR to speed up random sampling")
    } else {
      use.wrswoR <- FALSE
      message("### You are running advanced sampling. Installing the package \"wrswoR\" might result in much shorter runtimes.")
    }
  }

    # scale data if needed
  if(scale){
    message("### scaling input data...")
    x <- scale(x)
    # save the mean and stdev of the scaling
    x.scale.center <- attr(x = x, which = "scaled:center")
    x.scale.scale <- attr(x = x, which = "scaled:scale")
  }

  # make dir if needed
  if(!is.null(dir.randomization)){
    if(!dir.exists(dir.randomization))
      dir.create(dir.randomization)
  }



  ### get reference probabilities "Q"
  # using all points, set grid.points
  # from those, estimate Q
  # normalize to sum to 1
  message("### deciding grid points...")
  grid.coord <- get_grid_points(input=x, method=grid.method, grid.points=grid.points)

  # add another warning for the case that the number of grid.points was changed
  if(nrow(grid.coord) != grid.points){
    warning("The number of grid points was changed from ",grid.points," to ",nrow(grid.coord))
    grid.points <- nrow(grid.coord)
  }

  dist.to.grid <- get_dist_two_sets(x,grid.coord)

  # process the distances to a suitable density contribution
  # first, set bandwidth
  # bandwidth <- sqrt(sum((apply(x, 2, default_bandwidth.nrd)) ^ 2))
  # message("### using new bandwidth definition...")
  bandwidth <- median(apply(dist.to.grid,1,min))
  dist.to.grid.norm <- dist.to.grid / bandwidth
  density.contributions <-
    exp(-dist.to.grid.norm * dist.to.grid.norm / 2)

  if(is.null(use.advanced.sampling)){
    Q <- apply(density.contributions,2,sum)
  } else {
    Q <- apply(density.contributions*use.advanced.sampling,2,sum)
  }
  pseudo <- 1e-300 # quantile(Q[Q>0],0.01)
  Q <- Q + pseudo
  Q <- Q / sum(Q) # normalize to sum to 1



  ### get probabilities "P" for each group ("F" and "T")
  ### this has to be done for every gene X

  # for class "F" points, and for "T" points, separately
  # get 2D densities (using above grid points, limits, bandwidths)
  # add pseudocount to densities to avoid Inf problems
  # normalize to sum to 1
  # get D_KL (or relative entropy) of this P vs reference Q
  message("### calculating Kullback-Leibler divergences...")
  D_KL.observed <- rep(0,count.genes)
  class.types = c(FALSE, TRUE)
  pb <- txtProgressBar(min = 0, max = count.genes, style = 3, file = stderr()) # progress bar
  if(is.matrix(detection)){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL_highD(classes=detection[i,], density.contributions = density.contributions, reference.prob = Q, pseudo = pseudo)
      setTxtProgressBar(pb, i) # progress bar
    }
  } else if(inherits(detection, "lgRMatrix")){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL_highD(classes=extract_row_lgRMatrix(detection,i), density.contributions = density.contributions, reference.prob = Q, pseudo = pseudo)
      setTxtProgressBar(pb, i) # progress bar
    }
  }
  close(pb) # progress bar
  # return the sum of D_KL for "F" and "T"
  # store this value for each gene X



  ### use randomization to estimate expected values
  ### and p values for D_KL


  # Randomized D_KL values depend on the number of "F" and "T" cases
  # Therefore, for all observed numbers of "T" cases (from 0 to 100%):
  # - do x randomizations and get their D_KL
  # - get mean and SD per fraction,
  # - use those to estimate p values

  message("### performing randomizations...")
  T.counts <- Matrix::rowSums(detection)
  T.counts.unique <- sort(unique(T.counts))
  T.counts.unique.no <- length(T.counts.unique)
  p.vals <- rep(NA,nrow(detection))
  p.types <- rep(NA,nrow(detection))
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

  pb <- txtProgressBar(min = 0, max = T.counts.to.select, style = 3, file = stderr()) # progress bar

  # taking if - else statements outside of the for loop
  if(!is.null(use.advanced.sampling)){
    sampling.probs <- use.advanced.sampling

    for(i in 1:T.counts.to.select){
      setTxtProgressBar(pb, i) # progress bar

      T.count <- T.counts.selected[i]

      D_KL.randomized <- rep(NA,randomization.count)
      for(r in 1:randomization.count){
        # using sampling according to the number of genes expressed in each cell
        # pick cells according to the number of genes they express
        if(use.wrswoR){
          samp <- wrswoR::sample_int_expj(count.cells, size=T.count, prob=sampling.probs)
        } else {
          samp <- sample(x=count.cells, prob=sampling.probs, size=T.count, replace = FALSE)
        }
        # turn into T or F
        classes <- is.element(1:count.cells,samp)
        D_KL.randomized[r] <- get_D_KL_highD(classes=classes,
                                             density.contributions = density.contributions, reference.prob = Q, pseudo = pseudo)
      }
      all.D_KL.randomized[i,] <- D_KL.randomized
    }# end for all T counts to select
  } else {

    for(i in 1:T.counts.to.select){
      setTxtProgressBar(pb, i) # progress bar

      T.count <- T.counts.selected[i]

      D_KL.randomized <- rep(NA,randomization.count)
      vector.to.randomize <- c(rep(TRUE, T.count), rep(FALSE, ncol(detection)-T.count))
      for(r in 1:randomization.count){
        # using default sampling
        D_KL.randomized[r] <- get_D_KL_highD(classes=sample(x = vector.to.randomize),
                                             density.contributions = density.contributions, reference.prob = Q, pseudo = pseudo)

      }
      all.D_KL.randomized[i,] <- D_KL.randomized
    }# end for all T counts to select
  }# end if else
  close(pb) # progress bar

  message("### estimating p-values...")
  p.vals <- get_log_p_D_KL(T.counts = T.counts, D_KL.observed = D_KL.observed, D_KL.randomized = all.D_KL.randomized, output.dir = dir.randomization)

  # bonferroni correction for multiple testing
  p.adjs <- p.vals + log10(length(p.vals))
  p.adjs[p.adjs>0] <- 0 # p values should be at most 1; so log10 should be <= 0

  if(!is.null(dir.randomization)){
    message("### writing randomized Kullback-Leibler divergences to file...")
    outputfile.randomized.D_KL <- paste0(dir.randomization,"/random.D_KL.csv")
    write.csv(file=outputfile.randomized.D_KL,all.D_KL.randomized)
  }

  # prepare grid coordinates to return
  # if input data was scaled, the grid points have to be re-scaled
  # else nothing has to be done
  if(scale)
    grid.coord <- grid.coord*rep(x.scale.scale,each=nrow(grid.coord)) + rep(x.scale.center,each=nrow(grid.coord))

  message("### returning result...")
  # prepare the 'haystack' object to return
  res <- list(
    results = data.frame(
      D_KL = D_KL.observed,
      log.p.vals = p.vals,
      log.p.adj = p.adjs,
      T.counts = T.counts,
      row.names = row.names(detection)
    ),
    grid.coordinates = grid.coord
  )
  class(res) <- "haystack"
  res
}
