
#' The main Haystack function, for higher-dimensional spaces and continuous expression levels.
#'
#' @param x Coordinates of cells in a 2D or higher-dimensional space. Rows represent cells, columns the dimensions of the space.
#' @param expression a matrix with expression data of genes (rows) in cells (columns)
#' @param grid.points An integer specifying the number of centers (gridpoints) to be used for estimating the density distributions of cells. Default is set to 100.
#' @param weights.advanced.Q (Default: NULL) Optional weights of cells for calculating a weighted distribution of expression.
#' @param dir.randomization If NULL, no output is made about the random sampling step. If not NULL, files related to the randomizations are printed to this directory.
#' @param scale Logical (default=TRUE) indicating whether input coordinates in x should be scaled to mean 0 and standard deviation 1.
#' @param grid.method The method to decide grid points for estimating the density in the high-dimensional space. Should be "centroid" (default) or "seeding".
#' @param randomization.count Number of randomizations to use. Default: 100
#' @param n.genes.to.randomize Number of genes to use in randomizations. Default: 100
#'
#' @return An object of class "haystack", including the results of the analysis, and the coordinates of the grid points used to estimate densities.
#' @export
#'
#' @examples
#' # I need to add some examples.
#' # A toy example will be added too.
haystack_continuous_highD = function(x, expression, grid.points = 100, weights.advanced.Q=NULL, dir.randomization = NULL, scale=TRUE, grid.method="centroid", randomization.count = 100, n.genes.to.randomize = 100){
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
  if(!is.numeric(grid.points))
    stop("The value of 'grid.points' must be a numeric")
  if(grid.points >= ncol(expression))
    stop("The number of grid points appears to be very high (higher than the number of cells). You can set the number of grid points using the 'grid.points' parameter.")
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
  }# end for all T counts to select
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




#' The main Haystack function, for 2-dimensional spaces and continuous expression levels.
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param expression a matrix with expression data of genes (rows) in cells (columns)
#' @param weights.advanced.Q (Default: NULL) Optional weights of cells for calculating a weighted distribution of expression.
#' @param dir.randomization If NULL, no output is made about the random sampling step. If not NULL, files related to the randomizations are printed to this directory.
#' @param randomization.count Number of randomizations to use. Default: 100
#' @param n.genes.to.randomize Number of genes to use in randomizations. Default: 100
#'
#' @return An object of class "haystack", including the results of the analysis, and the coordinates of the grid points used to estimate densities.
#' @export
#'
#' @examples
#' # I need to add some examples.
#' # A toy example will be added too.
haystack_continuous_2D = function(x, y, expression, weights.advanced.Q = NULL, dir.randomization = NULL, randomization.count = 100, n.genes.to.randomize = 100){
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
  }# end for all T counts to select
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







#' Calculates the Kullback-Leibler divergence between distributions for the high-dimensional continuous version of haystack.
#'
#' @param weights A numerical vector with expression values of a gene.
#' @param density.contributions A matrix of density contributions of each cell (rows) to each center point (columns).
#' @param reference.prob A reference distribution to calculate the divergence against.
#' @param pseudo A pseudocount, used to avoid log(0) problems.
#'
#' @return A numerical value, the Kullback-Leibler divergence
get_D_KL_continuous_highD = function(weights, density.contributions, reference.prob, pseudo = 0){

  # the reference distribution Q of cells in the space
  Q <- reference.prob

  # calculating the Kullback-Leibler divergence of the distribution
  # of expression vs reference distribution Q
  P <- apply(density.contributions*weights, 2, sum)
  P <- P + pseudo
  P <- P / sum(P)
  D_KL <- sum(P * log(P/Q))

  D_KL
}


#' Calculates the Kullback-Leibler divergence between distributions for the 2-dimensional continuous version of haystack.
#'
#' @param weights A numerical vector with expression values of a gene.
#' @param parameters Parameters of the analysis, as set by function 'get_parameters_haystack'
#' @param reference.prob A reference distribution to calculate the divergence against.
#' @param pseudo A pseudocount, used to avoid log(0) problems.
#'
#' @return A numerical value, the Kullback-Leibler divergence
get_D_KL_continuous_2D = function(weights, parameters, reference.prob, pseudo){

  # the reference distribution Q of cells in the 2D plot
  Q <- reference.prob

  # calculating the Kullback-Leibler divergence of the distribution
  # of expression vs reference distribution Q
  density <- kde2d_faster(dens.x=parameters$dens.x*weights[col(parameters$dens.x)], # multiply cols of dens.x by weights
                          dens.y=parameters$dens.y*weights[col(parameters$dens.y)]) # same for dens.y

  # only decide pseudo if no value was given as input
  # (in total this takes some time)
  if(missing(pseudo))
    pseudo <- quantile(density[density>0],0.01)

  density <- density + pseudo
  P <- density / sum(density)

  D_KLs <- sum(P * log(P/Q))
  sum(D_KLs)
}

#' Estimates the significance of the observed Kullback-Leibler divergence by comparing to randomizations for the continuous version of haystack.
#'
#' @param D_KL.observed A vector of observed Kullback-Leibler divergences.
#' @param D_KL.randomized A matrix of Kullback-Leibler divergences of randomized datasets.
#' @param all.coeffVar Coefficients of variation of all genes. Used for fitting the Kullback-Leibler divergences.
#' @param train.coeffVar Coefficients of variation of genes that will be used for fitting the Kullback-Leibler divergences.
#' @param output.dir Optional parameter. Default is NULL. If not NULL, some files will be written to this directory.
#'
#' @return A vector of log10 p values, not corrected for multiple testing using the Bonferroni correction.
get_log_p_D_KL_continuous = function(D_KL.observed, D_KL.randomized, all.coeffVar, train.coeffVar, output.dir = NULL){

  # prepare data for fitting
  log2.D_KL <- log2(D_KL.randomized)
  D_KL.log.mean <- apply(log2.D_KL,1,mean)
  D_KL.log.sd   <- apply(log2.D_KL,1,sd)


  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  # first, for the mean D_KL
  x <- log(train.coeffVar)
  y <- D_KL.log.mean
  degree <- 3 # my simulations show that degree 3 works better than 1 or 2
  df <- 6 # my simulations show that this works best for n = 100. If n is higher, then df = 7 or 8 might become slightly better

  # set boundary knots to be slightly outside range of predictor values
  range.x <- range(x)
  r.x <- range.x[2]-range.x[1]
  boundary.knots <- c(range.x[1] - 0.01*r.x, range.x[2] + 0.01*r.x)

  model <- lm(y ~ bs(x, df=df,Boundary.knots=boundary.knots, degree=degree))
  #summary(model)

  if(!is.null(output.dir)){
    outfile <- paste0(output.dir,"/fit_logCoeffVar_vs_meanLogD_KL_df",df,"_degree",degree,".pdf")
    pdf(outfile)
    x.seq <- seq(range.x[1],range.x[2],length.out=1000)
    fitted.y <- predict(model, data.frame(x=x.seq), type="response")
    plot(x, y)
    points(x.seq,fitted.y,col="red", type="l")
    dev.off()
  }

  # fitted values for all cases
  fitted.D_KL_log.mean <- predict(model, data.frame(x=log(all.coeffVar)), type="response")


  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  # second, for the SD of the D_KL
  x <- train.coeffVar # not log in this case!
  y <- D_KL.log.sd
  degree <- 3 # my simulations show that degree 3 works better than 1 or 2
  df <- 5 # my simulations show that 5 or 6 works best for n = 100. If n is higher, then df = 7 or 8 might become slightly better

  # set boundary knots to be slightly outside range of predictor values
  range.x <- range(x)
  r.x <- range.x[2]-range.x[1]
  boundary.knots <- c(range.x[1] - 0.01*r.x, range.x[2] + 0.01*r.x)

  model <- lm(y ~ bs(x, df=df,Boundary.knots=boundary.knots, degree=degree))
  #summary(model)

  if(!is.null(output.dir)){
    outfile <- paste0(output.dir,"/fit_CoeffVar_vs_SdLogD_KL_df",df,"_degree",degree,".pdf")
    pdf(outfile)
    x.seq <- seq(range.x[1],range.x[2],length.out=1000)
    fitted.y <- predict(model, data.frame(x=x.seq), type="response")
    plot(x, y)
    points(x.seq,fitted.y,col="red", type="l")
    dev.off()
  }

  # fitted values for all cases
  fitted.D_KL_log.sd <- predict(model, data.frame(x=all.coeffVar), type="response")


  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  # estimate p values
  fitted.log.p.vals <- pnorm(log2(D_KL.observed), mean = fitted.D_KL_log.mean, sd = fitted.D_KL_log.sd, lower.tail = FALSE, log.p = T)/log(10)

  fitted.log.p.vals
}
