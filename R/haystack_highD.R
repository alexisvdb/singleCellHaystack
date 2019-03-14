# new function used for highD version
# get distances between all cells and all centers, in first PCs space
#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
get_euclidean_distance = function(x,y){
  sqrt(sum((x-y)^2))
}

# new function used for highD version
#' Title
#'
#' @param set1
#' @param set2
#'
#' @return
#' @export
#'
#' @examples
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

  class.types = c(F,T)

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
    P <- apply(density.contributions[cl.subset,],2,sum)
    P <- P/sum(P)
    D_KL <- sum(P * log(P/Q))
    D_KLs[c] <- D_KL
  }
  sum(D_KLs)
}



#' Title
#'
#' @param x
#' @param detection
#' @param centers
#' @param use.advanced.sampling
#' @param dir.randomization
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
haystack_highD = function(x, detection, centers = 50, use.advanced.sampling=NULL, dir.randomization = NULL, scale=TRUE){

  ##################################
  ##################################
  ### COMMENTED OUT FOR THE MOMENT
  ##################################
  ##################################
  # # check input
  # if(!is.numeric(x))
  #   stop("'x' must be a numeric vector")
  # if(!is.numeric(y))
  #   stop("'y' must be a numeric vector")
  # if(length(x) != length(y))
  #   stop("'x' and 'y' must have the same length")
  # if(ncol(detection) != length(x))
  #   stop("The number of columns in 'detection' must be the same as the length of 'x'")
  # if(!is.null(use.advanced.sampling)){
  #   if(!is.numeric(use.advanced.sampling))
  #     stop("'use.advanced.sampling' must either be NULL or a numeric vector")
  #   if(length(use.advanced.sampling) != length(x))
  #     stop("The length of 'use.advanced.sampling' must be the same as that of 'x'")
  # }
  #
  # # warn about unusal input sizes
  # if(length(x) < 50)
  #   warning("The number of cells seems very low (",length(x),"). Check your input.")
  # if(nrow(detection) < 100)
  #   warning("The number of genes seems very low (",nrow(detection),"). Check your input.")
  ##################################
  ##################################

  # scale data if needed
  if(scale){
    x <- scale(x)
  }

  # make dir if needed
  if(!is.null(dir.randomization)){
    if(!dir.exists(dir.randomization))
      dir.create(dir.randomization)
  }

  count.cells <- ncol(detection)
  count.genes <- nrow(detection)


  ### get reference probabilities "Q"

  # using all points, set centers
  # get distance between all points and all centers
  # set a reasonable bandwidth distance
  # get density contributions of each point to each center
  # from those, estimate Q
  # normalize to sum to 1

  center.n <- centers # 100 centers might be reasonable
  res.kmeans <- kmeans(x, centers=center.n, iter.max = 10, nstart = 10)
  center.coord <- res.kmeans$centers

  dist.to.centers <- get_dist_two_sets(x,center.coord)

  # process the distances to a suitable density contribution
  # first, set bandwidth
  bandwidth <- sqrt( sum( (apply(x,2,default_bandwidth.nrd)/4)^2 ) )
  dist.to.centers.norm <- dist.to.centers/bandwidth
  density.contributions <- exp(- dist.to.centers.norm*dist.to.centers.norm/2)

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
  message("### calculating Kulback-Leibler divergences...")
  D_KL.observed <- rep(0,count.genes)
  class.types = c(F,T)
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
        samp <- sample(x=count.cells, prob=sampling.probs, size=T.count, replace = F)
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
      vector.to.randomize <- c(rep(T,T.count),rep(F,ncol(detection)-T.count))
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
      row.names = row.names(detection)
    )
  )
  class(res) <- "haystack"
  res
}
