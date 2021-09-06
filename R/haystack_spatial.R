haystack_continuous_spatial = function(coordinates, expression, weights.advanced.Q=NULL, dir.randomization = NULL, scale=FALSE, randomization.count = 100, n.genes.to.randomize = 100){


  haystack_continuous_highD2(x = as.matrix(coordinates),
                             expression = expression,
                             grid.points = coordinates, # we are giving the coordinates here, not a number
                             weights.advanced.Q=weights.advanced.Q,
                             dir.randomization = dir.randomization,
                             scale=FALSE, # no need to scale, we are giving coordinates
                             grid.method="input",     # given as input
                             randomization.count = randomization.count,
                             n.genes.to.randomize = n.genes.to.randomize)

  x <- coordinates[,1]
  y <- coordinates[,2]
  haystack_continuous_2Dbis(x = x,
                            y = y,
                            expression = expression,
                            weights.advanced.Q=weights.advanced.Q,
                            dir.randomization = dir.randomization,
                            randomization.count = randomization.count,
                            n.genes.to.randomize = n.genes.to.randomize)

}



haystack_continuous_2Dbis = function(x, y, expression, weights.advanced.Q = NULL, dir.randomization = NULL, randomization.count = 100, n.genes.to.randomize = 100){
  message("### calling haystack_continuous_2D()...")
  message("### Using ",randomization.count," randomizations...")
  message("### Using ",n.genes.to.randomize," genes to randomize...")

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(!is.numeric(y))
    stop("'y' must be a numeric vector")
  if(length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if(!is.matrix(expression) && ! inherits(expression, "dgCMatrix") && ! inherits(expression, "dgRMatrix"))
    stop("'expression' must be a matrix, dgCMatrix, or dgRMatrix")
  if(ncol(expression) != length(x))
    stop("The number of columns in 'expression' must be the same as the length of 'x'")
  if(!is.null(weights.advanced.Q)){
    if(!is.numeric(weights.advanced.Q))
      stop("'weights.advanced.Q' must either be NULL or a numeric vector")
    if(length(weights.advanced.Q) != nrow(x))
      stop("The length of 'weights.advanced.Q' must be the same as the number of rows in 'x'")
  }

  # if expression is a dgCMatrix, convert it to a dgRMatrix
  if(inherits(expression, "dgCMatrix")){
    message("### converting expression data from dgCMatrix to dgRMatrix")
    expression <- as(expression, "RsparseMatrix")
  }

  count.cells <- ncol(expression)
  count.genes <- nrow(expression)
  if(n.genes.to.randomize > count.genes){
    warning("Number of genes to randomize (",n.genes.to.randomize,") is higher than the number of genes in the data.")
    warning("Seeting number of genes to randomize to ",count.genes)
    n.genes.to.randomize <- count.genes
  }

  # warn about unusual input sizes
  if(count.cells < 50)
    warning("The number of cells seems very low (",count.cells,"). Check your input.")
  if(count.cells > 10000)
    message("You are running haystack_2D on a large number of cells (",count.cells,"). Use method 'highD' for shorter runtimes.")
  if(count.genes < 100)
    warning("The number of genes seems very low (",count.genes,"). Check your input.")


  # make dir if needed
  if(!is.null(dir.randomization)){
    if(!dir.exists(dir.randomization))
      dir.create(dir.randomization)
  }

  ### get reference probabilities "Q"
  # using all points, get number of grid points, limits, and bandwidths
  # get 2D densities using all points (using above grid points, limits, bandwidths)
  # add pseudocount to densities to avoid Inf problems
  # normalize to sum to 1
  message("### setting parameters...")
  parameters <- get_parameters_haystack(x,y)

  # should be done using the default value for use.advanced.sampling = NULL
  # ref is a list that contains Q and the pseudo count that was used
  # ref$Q and red$pseudo
  ref <- get_reference(parameters, use.advanced.sampling = NULL)

  ### get density "P"
  ### this has to be one for every gene X
  ### now this is based on CONTINUOUS expression values
  ### and no longer on binary expression values (TRUE vs FALSE)


  # get densities (using above grid points, limits, bandwidths)
  # add pseudocount to densities to avoid Inf problems
  # normalize to sum to 1
  # get D_KL (or relative entropy) of this P vs reference Q
  message("### calculating Kullback-Leibler divergences...")
  D_KL.observed <- rep(0,count.genes)
  pb <- txtProgressBar(min = 0, max = count.genes, style = 3, file = stderr()) # progress bar
  if(is.matrix(expression)){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL_continuous_2D(weights=expression[i,], parameters=parameters, reference.prob=ref$Q, pseudo=ref$pseudo)
      setTxtProgressBar(pb, i) # progress bar
    }
  } else if(inherits(expression, "lgRMatrix")){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL_continuous_2D(weights=extract_row_lgRMatrix(expression,i), parameters=parameters, reference.prob=ref$Q, pseudo=ref$pseudo)
      setTxtProgressBar(pb, i) # progress bar
    }
  }
  close(pb) # progress bar


  ### use randomization to estimate expected values
  ### and p values for D_KL
  # for the mean D_KL:  the mean(log(D_KL)) can be modeled in function of log(coeffVar)
  # for the SD of D_KL: the sd(log(D_KL))   can be modeled in function of coeffVar
  # we need the Coefficient of Variation (coeffVar) ---> from the input data (expression)
  # we need the D_KL values of (a small set of) randomized genes ---> do randomizations

  # process expression data
  expr.mean <- apply(expression,1,mean) + 1e-300
  expr.sd   <- apply(expression,1,sd) + 1e-300
  coeffVar  <- expr.sd/expr.mean # coefficient of variation

  message("### performing randomizations...")

  o <- order(coeffVar)
  genes.to.randomize <- o[floor(seq(1,count.genes, length.out = n.genes.to.randomize))]

  # - do x randomizations and get their D_KL values
  # - get mean and SD of D_KL values,
  # - do fitting in function of coeffVar
  # - use those results to estimate p values of all genes

  # to store randomization results
  all.D_KL.randomized <- matrix(NA,nrow=n.genes.to.randomize, ncol=randomization.count)

  pb <- txtProgressBar(min = 0, max = n.genes.to.randomize, style = 3, file = stderr()) # progress bar

  for(i in 1:n.genes.to.randomize){
    setTxtProgressBar(pb, i) # progress bar

    D_KL.randomized <- rep(NA,randomization.count)
    vector.to.randomize <- expression[genes.to.randomize[i],]
    for(r in 1:randomization.count){
      # using default sampling
      D_KL.randomized[r] <- get_D_KL_continuous_2D(
        weights=sample(x=vector.to.randomize),
        parameters=parameters, reference.prob=ref$Q, pseudo=ref$pseudo
      )

    }
    all.D_KL.randomized[i,] <- D_KL.randomized
  }# end for all gene to randomize
  close(pb) # progress bar


  message("### estimating p-values...")
  p.vals <- get_log_p_D_KL_continuous(D_KL.observed = D_KL.observed,
                                      D_KL.randomized = all.D_KL.randomized,
                                      all.coeffVar = coeffVar,
                                      train.coeffVar = coeffVar[genes.to.randomize])

  # bonferroni correction for multiple testing
  p.adjs <- p.vals + log10(length(p.vals))
  p.adjs[p.adjs>0] <- 0 # p values should be at most 1; so log10 should be <= 0

  if(!is.null(dir.randomization)){
    message("### writing randomized Kullback-Leibler divergences to file...")
    outputfile.randomized.D_KL <- paste0(dir.randomization,"/random.D_KL.csv")
    write.csv(file=outputfile.randomized.D_KL,all.D_KL.randomized)
  }

  message("### returning result...")
  # prepare the 'haystack' object to return
  res <- list(
    results = data.frame(
      D_KL = D_KL.observed,
      log.p.vals = p.vals,
      log.p.adj = p.adjs,
      row.names = row.names(expression)
    )
  )
  class(res) <- "haystack"
  res
}







# a copy of haystack_continuous_highD on which I made smallish adjustments
haystack_continuous_highD2 = function(x, expression, grid.points = 100, weights.advanced.Q=NULL, dir.randomization = NULL, scale=TRUE, grid.method="centroid", randomization.count = 100, n.genes.to.randomize = 100){
  message("### calling haystack_continuous_highD()...")
  message("### Using ",randomization.count," randomizations...")
  message("### Using ",n.genes.to.randomize," genes to randomize...")

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric matrix")
  if(!is.matrix(x))
    stop("'x' must be a numeric matrix")
  if(ncol(x) < 2)
    stop("'x' must have at least 2 columns")
  if(!is.matrix(expression) && ! inherits(expression, "dgCMatrix") && ! inherits(expression, "dgRMatrix"))
    stop("'expression' must be a matrix, dgCMatrix, or dgRMatrix")
  if(ncol(expression) != nrow(x))
    stop("The number of columns in 'expression' must be the same as the rows in 'x'")
#  if(!is.numeric(grid.points))
#    stop("The value of 'grid.points' must be a numeric")
#  if(grid.points >= ncol(expression))
#    stop("The number of grid points appears to be very high (higher than the number of cells). You can set the number of grid points using the 'grid.points' parameter.")
  if(!is.null(weights.advanced.Q)){
    if(!is.numeric(weights.advanced.Q))
      stop("'weights.advanced.Q' must either be NULL or a numeric vector")
    if(length(weights.advanced.Q) != nrow(x))
      stop("The length of 'weights.advanced.Q' must be the same as the number of rows in 'x'")
  }

  # if expression is a dgCMatrix, convert it to a dgRMatrix
  if(inherits(expression, "dgCMatrix")){
    message("### converting expression data from dgCMatrix to dgRMatrix")
    expression <- as(expression, "RsparseMatrix")
  }

  count.cells <- ncol(expression)
  count.genes <- nrow(expression)
  if(n.genes.to.randomize > count.genes){
    warning("Number of genes to randomize (",n.genes.to.randomize,") is higher than the number of genes in the data.")
    warning("Seeting number of genes to randomize to ",count.genes)
    n.genes.to.randomize <- count.genes
  }

  # warn about unusal input sizes
  if(nrow(x) < 50)
    warning("The number of cells seems very low (",nrow(x),"). Check your input.")
  if(nrow(expression) < 100)
    warning("The number of genes seems very low (",nrow(expression),"). Check your input.")

  # warn about extreme values for 'grid.points'
  if(grid.points < 10)
    warning("The value of 'grid.points' appears to be very low (<10). You can set the number of grid points using the 'grid.points' parameter.")
  if(grid.points > count.cells/10)
    warning("The value of 'grid.points' appears to be very high (> No. of cells / 10). You can set the number of grid points using the 'grid.points' parameter.")

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
  if(grid.method == "input"){
    grid.coord <- x
  } els {
    grid.coord <- get_grid_points(input=x, method=grid.method, grid.points=grid.points)

    # add another warning for the case that the number of grid.points was changed
    if(nrow(grid.coord) != grid.points){
      warning("The number of grid points was changed from ",grid.points," to ",nrow(grid.coord))
      grid.points <- nrow(grid.coord)
    }

  }


  dist.to.grid <- get_dist_two_sets(x,grid.coord)

  # process the distances to a suitable density contribution
  # first, set bandwidth
  # bandwidth <- sqrt(sum((apply(x, 2, default_bandwidth.nrd)) ^ 2))
  # message("### using new bandwidth definition...")
  if(grid.method == "input"){
    tmp <- dist.to.grid
    tmp[tmp==0] <- NA
    bandwidth <- median(apply(tmp,1,function(x) min(x, na.rm = TRUE)))
  } else {
    bandwidth <- median(apply(dist.to.grid,1,min))
  }
  dist.to.grid.norm <- dist.to.grid / bandwidth
  density.contributions <-
    exp(-dist.to.grid.norm * dist.to.grid.norm / 2)

  if(is.null(weights.advanced.Q)){
    Q <- apply(density.contributions,2,sum)
  } else {
    Q <- apply(density.contributions*weights.advanced.Q,2,sum)
  }
  pseudo <- 1e-300 # quantile(Q[Q>0],0.01)
  Q <- Q + pseudo
  Q <- Q / sum(Q) # normalize to sum to 1



  ### get density "P"
  ### now this is based on CONTINUOUS expression values
  ### and no longer on binary expression values (TRUE vs FALSE)

  # get densities (using above grid points, limits, bandwidths)
  # add pseudocount to densities to avoid Inf problems
  # normalize to sum to 1
  # get D_KL (or relative entropy) of this P vs reference Q
  message("### calculating Kullback-Leibler divergences...")
  D_KL.observed <- rep(0,count.genes)
  pb <- txtProgressBar(min = 0, max = count.genes, style = 3, file = stderr()) # progress bar
  if(is.matrix(expression)){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL_continuous_highD(weights=expression[i,], density.contributions = density.contributions, reference.prob = Q, pseudo = pseudo)
      setTxtProgressBar(pb, i) # progress bar
    }
  } else if(inherits(expression, "dgRMatrix")){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL_continuous_highD(weights=extract_row_dgRMatrix(expression,i), density.contributions = density.contributions, reference.prob = Q, pseudo = pseudo)
      setTxtProgressBar(pb, i) # progress bar
    }
  }
  close(pb) # progress bar


  ### use randomization to estimate expected values
  ### and p values for D_KL
  # for the mean D_KL:  the mean(log(D_KL)) can be modeled in function of log(coeffVar)
  # for the SD of D_KL: the sd(log(D_KL))   can be modeled in function of coeffVar
  # we need the Coefficient of Variation (coeffVar) ---> from the input data (expression)
  # we need the D_KL values of (a small set of) randomized genes ---> do randomizations

  # process expression data
  expr.mean <- apply(expression,1,mean) + 1e-300
  expr.sd   <- apply(expression,1,sd) + 1e-300
  coeffVar  <- expr.sd/expr.mean # coefficient of variation

  message("### performing randomizations...")

  o <- order(coeffVar)
  genes.to.randomize <- o[floor(seq(1,count.genes, length.out = n.genes.to.randomize))]

  # - do x randomizations and get their D_KL values
  # - get mean and SD of D_KL values,
  # - do fitting in function of coeffVar
  # - use those results to estimate p values of all genes

  # to store randomization results
  all.D_KL.randomized <- matrix(NA,nrow=n.genes.to.randomize, ncol=randomization.count)

  pb <- txtProgressBar(min = 0, max = n.genes.to.randomize, style = 3, file = stderr()) # progress bar

  for(i in 1:n.genes.to.randomize){
    setTxtProgressBar(pb, i) # progress bar

    D_KL.randomized <- rep(NA,randomization.count)
    vector.to.randomize <- expression[genes.to.randomize[i],]
    for(r in 1:randomization.count){
      # using default sampling
      D_KL.randomized[r] <- get_D_KL_continuous_highD(
        weights=sample(x=vector.to.randomize),
        density.contributions = density.contributions, reference.prob = Q, pseudo = pseudo
      )

    }
    all.D_KL.randomized[i,] <- D_KL.randomized
  }# end for all gene to randomize
  close(pb) # progress bar


  message("### estimating p-values...")
  p.vals <- get_log_p_D_KL_continuous(D_KL.observed = D_KL.observed,
                                      D_KL.randomized = all.D_KL.randomized,
                                      all.coeffVar = coeffVar,
                                      train.coeffVar = coeffVar[genes.to.randomize])

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
      row.names = row.names(expression)
    ),
    grid.coordinates = grid.coord #,
    #all.D_KL.randomized = all.D_KL.randomized
  )
  class(res) <- "haystack"
  res
}

