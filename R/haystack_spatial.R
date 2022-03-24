#' The main Haystack function for spatial transcriptomics data and continuous expression levels.
#'
#' @param coordinates Coordinates of cells in spatial transcriptomics data.
#' @param expression a matrix with expression data of genes (rows) in cells or spots (columns)
#' @param weights.advanced.Q (Default: NULL) Optional weights of cells for calculating a weighted distribution of expression.
#' @param dir.randomization If NULL, no output is made about the random sampling step. If not NULL, files related to the randomizations are printed to this directory.
#' @param scale Logical (default=FALSE) indicating whether input coordinates should be scaled to mean 0 and standard deviation 1.
#' @param randomization.count Number of randomizations to use. Default: 100
#' @param n.genes.to.randomize Number of genes to use in randomizations. Default: 100
#'
#' @return An object of class "haystack", including the results of the analysis.
#' @export
#'
#' @examples
#' # I need to add some examples.
#' # A toy example will be added too.
haystack_continuous_spatial = function(coordinates, expression, weights.advanced.Q=NULL, dir.randomization = NULL, scale=FALSE, randomization.count = 100, n.genes.to.randomize = 100){

  message("### calling haystack_continuous_spatial2D()...")
  message("### Using ",randomization.count," randomizations...")
  message("### Using ",n.genes.to.randomize," genes to randomize...")

  # check input
  if(ncol(coordinates)>3){
    warning("Input coordinates appear to have ",ncol(coordinates)," dimensions. Spatial transcriptomics data typically has 2 or at most 3.")
    warning("Please make sure that running the haystack function for spatial transcriptomics data makes sense for your data.")
  }
  if(!is.numeric(coordinates) && !all(apply(coordinates, 2, is.numeric)))
    stop("'coordinates' must be a numeric matrix or data.frame")
  if(!is.matrix(coordinates) && !is.data.frame(coordinates))
    stop("'coordinates' must be a numeric matrix or data.frame")
  if(!is.matrix(expression) && ! inherits(expression, "dgCMatrix") && ! inherits(expression, "dgRMatrix"))
    stop("'expression' must be a matrix, dgCMatrix, or dgRMatrix")
  if(ncol(expression) != nrow(coordinates))
    stop("The number of columns in 'expression' must be the same as the number of rows in 'coordinates'")
  if(!is.null(weights.advanced.Q)){
    if(!is.numeric(weights.advanced.Q))
      stop("'weights.advanced.Q' must either be NULL or a numeric vector")
    if(length(weights.advanced.Q) != nrow(coordinates))
      stop("The length of 'weights.advanced.Q' must be the same as the number of rows in 'coordinates'")
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
    warning("The number of cells or spots seems very low (",count.cells,"). Check your input.")
  if(count.genes < 100)
    warning("The number of genes seems very low (",count.genes,"). Check your input.")

  # make dir if needed
  if(!is.null(dir.randomization)){
    if(!dir.exists(dir.randomization))
      dir.create(dir.randomization)
  }

  ### get reference probabilities "Q"
  # using all points also as grid points
  # get pairwise distances, bandwidths, and density contributions
  # based on those, get "Q"
  # add pseudocount to densities to avoid Inf problems
  # normalize to sum to 1

  message("### getting pairwise density contributions...")
  # I found that putting this in a function is faster than putting it here
  pairwise.densities <- get_parameters_haystack_spatial(coordinates = coordinates, bandwidth.multiplier = 1)

  # should be done using the default value for use.advanced.sampling = NULL
  # ref is a list that contains Q and the pseudo count that was used
  # ref$Q and red$pseudo
  ref <- get_reference_spatial(pairwise.densities, use.advanced.sampling = NULL)




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
      D_KL.observed[i] <- get_D_KL_continuous_spatial(weights=expression[i,], density.contributions=pairwise.densities, reference.prob=ref$Q, pseudo=0)
      setTxtProgressBar(pb, i) # progress bar
    }
  } else if(inherits(expression, "dgRMatrix")){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL_continuous_spatial(weights=extract_row_dgRMatrix(expression,i), density.contributions=pairwise.densities, reference.prob=ref$Q, pseudo=0)
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

    if(is.matrix(expression)){
      vector.to.randomize <- expression[genes.to.randomize[i],]
    } else if(inherits(expression, "dgRMatrix")){
      vector.to.randomize <- extract_row_dgRMatrix(expression,genes.to.randomize[i])
    }

    for(r in 1:randomization.count){
      # using default sampling
      D_KL.randomized[r] <- get_D_KL_continuous_spatial(
        weights=sample(x=vector.to.randomize),
        density.contributions = pairwise.densities, reference.prob=ref$Q, pseudo=0
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
  info <- p.vals$info
  info$mean$observed$feature <- rownames(expression[genes.to.randomize, ])
  info$sd$observed$feature <- rownames(expression[genes.to.randomize, ])
  p.vals <- p.vals$fitted

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
    ),
    info = list(
      method="continuous_spatial",
      randomization = info
    )
  )
  class(res) <- "haystack"
  res

}

#' Decides the pairwise density contributions between cells or spots, which will be used during the "Haystack" analysis.
#'
#' @param coordinates Coordinates of cells in spatial transcriptomics data.
#' @param bandwidth.multiplier A parameter to adjust the bandwitdh used for estimating the density contributions. Default = 1.
#'
#' @return Pairwise density contributions between all pairs of cell or spots.
get_parameters_haystack_spatial = function(coordinates, bandwidth.multiplier = 1){

  # get pairwise distances between spots
  pairwise.dists <- as.matrix(dist(coordinates))

  # set bandwidth
  tmp <- pairwise.dists
  tmp[tmp == 0] <- NA
  min.dists <- apply(tmp, 1, function(x) min(x, na.rm = TRUE))
  bandwidth <- median(min.dists) * bandwidth.multiplier

  # getting the distances (in units of bandwidths) between all cells and all grid points
  pairwise.dists.norm <- pairwise.dists / min.dists

  # densities based on a simplified Gaussian
  pairwise.densities <- exp(-0.5 * pairwise.dists.norm*pairwise.dists.norm)

  pairwise.densities
}




#' Get reference distribution for spatial transcriptomics data
#'
#' @param pairwise.densities Pairwise density contributions between all pairs of grid points
#' @param use.advanced.sampling If NULL naive sampling is used. If a vector is given (of length = no. of cells) sampling is done according to the values in the vector.
#'
#' @return A list with two components, Q for the reference distribution and pseudocounts pseudo.
#'
get_reference_spatial <- function(pairwise.densities, use.advanced.sampling = NULL) {
  density <- apply(pairwise.densities,1,sum)

  if (!is.null(use.advanced.sampling)) {
    density <- density * use.advanced.sampling
  }
  Q.tmp = density / sum(density)

  # there are no non-zero values now.
  # most spots are surrounded by other spots.
  # in the other functions, the pseudocount was used to avoid 0 values
  # here, I am not yet sure how to pick a proper psuedo count
  # in fact, I am not sure if a pseudocount is needed
  pseudo <- 0

  density <- density + pseudo

  Q = density / sum(density)
  list(Q = Q, pseudo = pseudo)
}


#' Calculates the Kullback-Leibler divergence between distributions for the spatial transcriptomics version of haystack.
#'
#' @param weights A numerical vector with expression values of a gene.
#' @param density.contributions A matrix of pairwise density contributions.
#' @param reference.prob A reference distribution to calculate the divergence against.
#' @param pseudo A pseudocount, used to avoid log(0) problems.
#'
#' @return A numerical value, the Kullback-Leibler divergence
get_D_KL_continuous_spatial = function(weights, density.contributions, reference.prob, pseudo = 0){

  # the reference distribution Q of cells in the space
  Q <- reference.prob

  # calculating the Kullback-Leibler divergence of the distribution
  # of expression vs reference distribution Q
  # P <- apply(density.contributions*weights, 2, sum)
  P <- colSums(density.contributions*weights)
  P <- P / sum(P)
  D_KL <- sum(P * log(P/Q))

  D_KL
}

