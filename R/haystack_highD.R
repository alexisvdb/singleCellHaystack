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
#' @param pseudo A pseudocount, used to avoid log(0) problems. Should not be needed, so default is 0.
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
    if(sum(cl.subset)==0)
      next
    P <- apply(density.contributions[cl.subset, , drop = FALSE], 2, sum)
    P <- P/sum(P)
    D_KL <- sum(P * log(P/Q))
    D_KLs[c] <- D_KL
  }
  sum(D_KLs)
}

#' Function for finding for a point in space the distance to the nearest point with higher density. See
#'
#' @param entries_n The number of points.
#' @param dist A distance object with the distances between all pairs oof points.
#' @param degree_counts A vector with the degree or density of each point.
#'
#' @return A vector with distances.
get.distance.to.nearest.higher.degree = function(entries_n=length(degree_counts), dist, degree_counts){

  # run through the distances between pairs of points
  # each time, check which of the pair has the highest degree
  # keep track of "closest" entries with higher degree
  dist_vector <- as.vector(dist)
  dist_n <- entries_n * (entries_n-1) / 2
  index_1 <- 1
  index_2 <- 2
  max_distance <- max(dist)
  distance_to_nearest_higher_degree <- apply(as.matrix(dist),1, max)
  for(i in 1:dist_n){
    c <- dist_vector[i]

    degree_1 <- degree_counts[index_1]
    degree_2 <- degree_counts[index_2]

    if(degree_1 > degree_2){
      if(c < distance_to_nearest_higher_degree[index_2]){
        distance_to_nearest_higher_degree[index_2] <- c
      }
    } else if (degree_2 > degree_1){
      if(c < distance_to_nearest_higher_degree[index_1]){
        distance_to_nearest_higher_degree[index_1] <- c
      }
    }

    # update indices
    index_2 <- index_2 + 1
    if(index_2 > entries_n){
      index_1 <- index_1 + 1
      index_2 <- index_1 + 1
    }

  }

  # return result
  distance_to_nearest_higher_degree

}# end function get.distance.to.nearest.higher.degree



#' A function to decide grid points in a higher-dimensional space
#'
#' @param input A numerical matrix with higher-dimensional coordinates (columns) of points (rows)
#' @param method The method to decide grid points. Should be "grid" (default) or "kmeans".
#' @param grid.points The number of grid points to return. Default is 50.
#'
#' @return Coordinates of grid points in the higher-dimensonal space.
get_grid_points = function(input, method="grid", grid.points = 50){

  if(method=="kmeans"){
    # perform k-means clustering and get the centers
    res.kmeans <- kmeans(input, centers=grid.points, iter.max = 10, nstart = 10)
    grid.coord <- res.kmeans$centers

  } else if(method=="grid"){
    # first decide a grid in all dimensions
    # then, "round" cells to the nearest grid point
    # from all gird points, keep those that are closest to many cells, and far away from other points

    input.dim <- ncol(input)

    # get the 10% and 90% points in each dimension
    input10 <- apply(input, 2, function(x) as.numeric(quantile(x,0.10)))
    input90 <- apply(input, 2, function(x) as.numeric(quantile(x,0.90)))
    # set bin size of each dimension
    grid.points.10.90 <- round(log(x=2000,base=input.dim)) # naively, this should result in about 2000 initial grid points in total
    input.bin.size <- (input90-input10)/grid.points.10.90

    # round to nearest grid point
    input.grid <- t(apply(input,1, function(a) round(a / input.bin.size) * input.bin.size))
    ### no idea why this transposes this!! So, for the moment I added t()

    # filter out grid points assigned to only 1 cell
    count.threshold.grid <- 1
    input.table <- data.table::as.data.table(input.grid)
    N <- nrow(input.table)
    input.table.count <- input.table[, .N, by = c(colnames(input))]
    input.table.subset <- subset(input.table.count, N>count.threshold.grid)

    degrees <- input.table.subset$N
    input.grid.candidates <- as.matrix(input.table.subset)[,-ncol(input.table.subset)]
    dist.input.grid.candidats <- dist(input.grid.candidates)

    dist.nearest.higher.degree <- get.distance.to.nearest.higher.degree(entries_n=nrow(input.grid.candidates), dist=dist.input.grid.candidats, degree_counts = degrees)
    # gamma <- degrees*dist.nearest.higher.degree # original version
    gamma <- log2(degrees)*dist.nearest.higher.degree # more stress on spreading out points
    # plot(sort(gamma, decreasing = T))

    center.indices <- order(gamma,decreasing = T)[1:grid.points]
    grid.coord <- input.grid.candidates[center.indices,]
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
#' @param grid.points An integer specifying the number of centers (gridpoints) to be used for estimating the density distributions of cells. Default is set to 50.
#' @param grid.method The method to decide grid points for estimating the density in the high-dimensional space. Should be "grid" (default) or "kmeans".
#'
#' @return An object of class "haystack", including the results of the analysis, and the coordinates of the grid points used to estimate densities.
#' @export
#'
#' @examples
#' # I need to add some examples.
#' # A toy example will be added too.
haystack_highD = function(x, detection, grid.points = 50, use.advanced.sampling=NULL, dir.randomization = NULL, scale=TRUE, grid.method="grid"){
  message("### calling haystack_highD()...")

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric matrix")
  if(!is.matrix(x))
    stop("'x' must be a numeric matrix")
  if(ncol(x) < 2)
    stop("'x' must have at least 2 columns")
  if(!is.matrix(detection))
    stop("'detection' must be a matrix")
  if(ncol(detection) != nrow(x))
    stop("The number of columns in 'detection' must be the same as the rows in 'x'")
  if(!is.numeric(grid.points))
    stop("The value of 'grid.points' must be a numeric")
  if(grid.points >= ncol(detection))
    stop("The number of grid.points appears to be very high; higher than the number of cells. Please check your input.")
  if(!is.null(use.advanced.sampling)){
    if(!is.numeric(use.advanced.sampling))
      stop("'use.advanced.sampling' must either be NULL or a numeric vector")
    if(length(use.advanced.sampling) != nrow(x))
      stop("The length of 'use.advanced.sampling' must be the same as the number of rows in 'x'")
  }
  if(!is.logical(scale) | length(scale) > 1)
    stop("The value of 'scale' must be either TRUE or FALSE")
  if(grid.method!="grid" & grid.method!="kmeans")
    stop("The value of 'grid.method' must be either 'grid' or 'kmeans'")

  count.cells <- ncol(detection)
  count.genes <- nrow(detection)

  # warn about unusal input sizes
  if(nrow(x) < 50)
    warning("The number of cells seems very low (",nrow(x),"). Check your input.")
  if(nrow(detection) < 100)
    warning("The number of genes seems very low (",nrow(detection),"). Check your input.")

  # warn about extreme values for 'grid.points'
  if(grid.points < 10)
    warning("The value of 'grid.points' appears to be very low (<10). Please check your input.")
  if(grid.points > count.cells/10)
    warning("The value of grid.points appears to be very high (> No. of cells / 10). Please check your input.")


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

  dist.to.grid <- get_dist_two_sets(x,grid.coord)

  # process the distances to a suitable density contribution
  # first, set bandwidth
  bandwidth <- sqrt( sum( (apply(x,2,default_bandwidth.nrd)/4)^2 ) )
  dist.to.grid.norm <- dist.to.grid/bandwidth
  density.contributions <- exp(- dist.to.grid.norm*dist.to.grid.norm/2)

  if(is.null(use.advanced.sampling)){
    Q <- apply(density.contributions,2,sum)
  } else {
    Q <- apply(density.contributions*use.advanced.sampling,2,sum)
  }
  # NO NEED FOR PSEUDOCOUNT HERE??
  Q <- Q/sum(Q) # normalize to sum to 1



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
  for(i in 1:count.genes){
    D_KL.observed[i] <- get_D_KL_highD(classes=detection[i,], density.contributions = density.contributions, reference.prob = Q)

    if(i%%1000==0)
      message(paste0("### ... ",i," values out of ",count.genes," done"))
  }
  # return the sum of D_KL for "F" and "T"
  # store this value for each gene X



  ### use randomization to estimate expected values
  ### and p values for D_KL


  # Randomized D_KL values depend on the number of "F" and "T" cases
  # Thereofre, for all observed numbers of "T" cases (from 0 to 100%):
  # - do x randomizations and get teir D_KL
  # - get mean and SD per fraction,
  # - use those to estimate p values

  message("### starting randomizations...")
  T.counts <- apply(detection,1,sum)
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
        samp <- sample(x=count.cells, prob=sampling.probs, size=T.count, replace = FALSE)
        # turn into T or F
        classes <- is.element(1:count.cells,samp)
        D_KL.randomized[r] <- get_D_KL_highD(classes=classes,
                                             density.contributions = density.contributions, reference.prob = Q)
      }
      all.D_KL.randomized[i,] <- D_KL.randomized
    }# end for all T counts to select
  } else {

    for(i in 1:T.counts.to.select){
      if(i%%10==0)
        message(paste0("### ... ",i," values out of ",T.counts.to.select," done"))

      T.count <- T.counts.selected[i]

      D_KL.randomized <- rep(NA,randomization.count)
      vector.to.randomize <- c(rep(TRUE, T.count), rep(FALSE, ncol(detection)-T.count))
      for(r in 1:randomization.count){
        # using default sampling
        D_KL.randomized[r] <- get_D_KL_highD(classes=sample(x = vector.to.randomize),
                                             density.contributions = density.contributions, reference.prob = Q)

      }
      all.D_KL.randomized[i,] <- D_KL.randomized
    }# end for all T counts to select
  }# end if else

  message("### estimating p-values...")
  p.vals <- get_log_p_D_KL(T.counts = T.counts, D_KL.observed = D_KL.observed, D_KL.randomized = all.D_KL.randomized, output.dir = dir.randomization)

  if(!is.null(dir.randomization)){
    message("### writing randomized Kullback-Leibler divergences to file...")
    outputfile.randomized.D_KL <- paste0(dir.randomization,"/random.D_KL.csv")
    write.csv(file=outputfile.randomized.D_KL,all.D_KL.randomized)
  }

  # prepare grid coordinates to return
  # if input data was scaled, the grid points have to be re-scaled
  # else nothing has to be done
  if(scale)
    grid.coord <- grid.coord*x.scale.scale + x.scale.center

  message("### returning result...")
  # prepare the 'haystack' object to return
  res <- list(
    results = data.frame(
      D_KL = D_KL.observed,
      log.p.vals = p.vals,
      T.counts = T.counts,
      row.names = row.names(detection)
    ),
    grid.coordinates = grid.coord
  )
  class(res) <- "haystack"
  res
}
