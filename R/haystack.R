#' Based on the MASS kde2d() function, but heavily simplified; it's just tcrossprod() now.
#'
#' @param dens.x Contribution of all cells to densities of the x-axis grid points.
#' @param dens.y Contribution of all cells to densities of the y-axis grid points.
#'
kde2d_faster = function (dens.x, dens.y){
  tcrossprod(dens.x, dens.y)
}

#' Default function given by function bandwidth.nrd in MASS. No changes were made to this function.
#'
#' @param x A numeric vector
#'
#' @return A suitable bandwith.
default_bandwidth.nrd = function(x){
  r <- quantile(x, c(0.25, 0.75))
  h <- (r[2] - r[1]) / 1.34
  1.06 * min(sqrt(var(x)), h) * length(x) ^ (-1 / 5)
}

#' Function that decides most of the parameters that will be during the "Haystack" analysis.
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param high.resolution Logical: should high resolution be used? Default is FALSE.
#'
#' @return A list containing various parameters to use in the analysis.
get_parameters_haystack = function(x, y, high.resolution = FALSE){

  # the bandwidths used for the kernel function
  bandwidth.x <- default_bandwidth.nrd(x)
  bandwidth.y <- default_bandwidth.nrd(y)
  bandwidths <- c(bandwidth.x,bandwidth.y)
  bandwidths <- bandwidths

  # I want to have a fixed number of bins between the 10th and 90th percentile in each dimension
  # this is to avoid bins getting squished because of a few outliers
  # the default number of bins between the 10% and 90% points in each dimension
  if(high.resolution) {
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
  # with an added padding the size of a bandwidth
  x.min <- min(x) - bandwidths[1]
  x.max <- max(x) + bandwidths[1]
  y.min <- min(y) - bandwidths[2]
  y.max <- max(y) + bandwidths[2]

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

  # getting the distances (in units of bandwidths) between all cells and all grid points
  gx <- seq.int(x.lim.min, x.lim.max, length.out = x.total.grid.points)
  gy <- seq.int(y.lim.min, y.lim.max, length.out = y.total.grid.points)
  ax <- outer(gx, x, "-")/bandwidths[1L]
  ay <- outer(gy, y, "-")/bandwidths[2L]

  # densities based on a simplified Gaussian
  dens.x <- matrix(exp(-0.5 * ax * ax), ncol = length(x))
  dens.y <- matrix(exp(-0.5 * ay * ay), ncol = length(x)) # length(x) should be the same as length(y), of course

  # return results
  list(
    grid.points = grid.points,
    limits      = limits,
    bandwidths  = bandwidths,
    dens.x      = dens.x,
    dens.y      = dens.y
  )
}

#' Calculates the Kullback-Leibler divergence between distributions.
#'
#' @param classes A logical vector. Values are T is the gene is expressed in a cell, F is not.
#' @param parameters Parameters of the analysis, as set by function 'get_parameters_haystack'
#' @param reference.prob A reference distribution to calculate the divergence against.
#' @param pseudo A pseudocount, used to avoid log(0) problems.
#'
#' @return A numerical value, the Kullback-Leibler divergence
get_D_KL = function(classes, parameters, reference.prob, pseudo){

  class.types = c(FALSE, TRUE)

  # the reference distribution Q of cells in the 2D plot
  Q <- reference.prob

  # calculating the Kullback-Leibler divergence of the distribution
  # of cells expressin and not expressing gene X vs reference distribution Q
  D_KLs <- rep(NA,length(class.types))
  for(c in 1:length(class.types)){
    cl <- class.types[c]
    cl.subset <- classes==cl
    if(sum(cl.subset)==0){ # this means a gene is expressed or not expressed in all cells; return 0
      #D_KLs[c] <- 0
      return(0)
    } else {
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
  }
  sum(D_KLs)
}

#' Estimates the significance of the observed Kullback-Leibler divergence by comparig to randomizations.
#'
#' @param T.counts The number of cells in which a gene is detected.
#' @param D_KL.observed A vector of observed Kullback-Leibler divergences.
#' @param D_KL.randomized A matrix of Kullback-Leibler divergences of randomized datasets.
#' @param output.dir Optional parameter. Default is NULL. If not NULL, some files will be written to this directory.
#'
#' @return A vector of log10 p values, not corrected for multiple testing using the Bonferroni correction.
get_log_p_D_KL = function(T.counts, D_KL.observed, D_KL.randomized, output.dir = NULL){

  t.points <- as.numeric(row.names(D_KL.randomized))
  dat.mean.log2 <- apply(log2(D_KL.randomized),1,mean)
  dat.sd.log2 <- apply(log2(D_KL.randomized),1,sd)

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

  # apply on actual observed values

  fitted.mean.log2 <- predict(model.mean.log2, data.frame(t.points=T.counts), type="response")
  fitted.sd.log2 <- predict(model.sd.log2, data.frame(t.points=T.counts), type="response")

  fitted.log.p.vals <- pnorm(log2(D_KL.observed), mean = fitted.mean.log2, sd = fitted.sd.log2, lower.tail = FALSE, log.p = T)/log(10)

  fitted.log.p.vals
}


#' Get reference distribution
#'
#' @param param Parameters of the analysis, as set by function 'get_parameters_haystack'
#' @param use.advanced.sampling If NULL naive sampling is used. If a vector is given (of length = no. of cells) sampling is done according to the values in the vector.
#'
#' @return A list with two components, Q for the reference distribution and pseudo.
#'
get_reference <- function(param, use.advanced.sampling = NULL) {
  if (is.null(use.advanced.sampling)) {
    density <- kde2d_faster(param$dens.x, param$dens.y)
    pseudo <- quantile(density[density > 0], 0.01)
  } else {
    density <- kde2d_faster(t(t(param$dens.x) * use.advanced.sampling),
                            t(t(param$dens.y)))
    density <- density / sum(density)
    pseudo <- quantile(density[density>0],0.01)
  }
  density <- density + pseudo

  Q = density / sum(density)
  list(Q = Q, pseudo = pseudo)
}

#' The main Haystack function, for 2-dimensional spaces.
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param detection A logical matrix showing which genes (rows) are detected in which cells (columns)
#' @param use.advanced.sampling If NULL naive sampling is used. If a vector is given (of length = no. of cells) sampling is done according to the values in the vector.
#' @param dir.randomization If NULL, no output is made about the random sampling step. If not NULL, files related to the randomizations are printed to this directory.
#'
#' @return An object of class "haystack"
#' @export
#'
#' @examples
#' # using the toy example of the singleCellHaystack package
#' # define a logical matrix with detection of each gene (rows) in each cell (columns)
#' dat.detection <- dat.expression > 1
#'
#' # running haystack in default mode
#' res <- haystack(dat.tsne, detection=dat.detection, method = "2D")
#' # list top 10 biased genes
#' show_result_haystack(res, n =10)
haystack_2D = function(x, y, detection, use.advanced.sampling=NULL, dir.randomization = NULL){
  message("### calling haystack_2D()...")

  # check input
  if(!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(!is.numeric(y))
    stop("'y' must be a numeric vector")
  if(length(x) != length(y))
    stop("'x' and 'y' must have the same length")
  if(!is.matrix(detection) && ! inherits(detection, "lgCMatrix") && ! inherits(detection, "lgRMatrix"))
    stop("'detection' must be a matrix, lgCMatrix, or lgRMatrix")
  if(ncol(detection) != length(x))
    stop("The number of columns in 'detection' must be the same as the length of 'x'")
  if(!is.null(use.advanced.sampling)){
    if(!is.numeric(use.advanced.sampling))
      stop("'use.advanced.sampling' must either be NULL or a numeric vector")
    if(length(use.advanced.sampling) != length(x))
      stop("The length of 'use.advanced.sampling' must be the same as that of 'x'")
  }

  count.cells <- ncol(detection)
  count.genes <- nrow(detection)

  # warn about unusual input sizes
  if(count.cells < 50)
    warning("The number of cells seems very low (",count.cells,"). Check your input.")
  if(count.cells > 10000)
    message("You are running haystack_2D on a large number of cells (",count.cells,"). Use method 'highD' for shorter runtimes.")
  if(count.genes < 100)
    warning("The number of genes seems very low (",count.genes,"). Check your input.")
  # if detection is a lgCMatrix, convert it to a lgRMatrix
  if(inherits(detection, "lgCMatrix")){
    message("### converting detection data from lgCMatrix to lgRMatrix")
    detection <- as(detection, "RsparseMatrix")
  }

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

  # if(is.null(use.advanced.sampling)){
  #   density <- kde2d_faster(dens.x=parameters$dens.x,dens.y=parameters$dens.y)
  #   pseudo <- quantile(density[density>0],0.01)
  # } else {
  #   density <- kde2d_faster(dens.x=t(t(parameters$dens.x)*use.advanced.sampling),
  #                           dens.y=t(t(parameters$dens.y)))
  #   density <- density / sum(density)
  #   pseudo <- quantile(density[density>0],0.01)
  # }
  # density <- density + pseudo
  # # heatmap(density, Rowv=NA,Colv=NA,scale="none")
  # Q <- density / sum(density)
  ref <- get_reference(parameters, use.advanced.sampling)

  ### get probabilities "P" for each group ("F" and "T")
  ### this has to be one for every gene X


  # for class "F" points, and for "T" points, separately
  # get 2D densities (using above grid points, limits, bandwidths)
  # add pseudocount to densities to avoid Inf problems
  # normalize to sum to 1
  # get D_KL (or relative entropy) of this P vs reference Q
  message("### calculating Kullback-Leibler divergences...")
  D_KL.observed <- rep(0,count.genes)
  pb <- txtProgressBar(min = 0, max = count.genes, style = 3, file = stderr()) # progress bar
  if(is.matrix(detection)){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL(classes=detection[i,], parameters=parameters, reference.prob=ref$Q, pseudo=ref$pseudo)
      setTxtProgressBar(pb, i) # progress bar
    }
  } else if(inherits(detection, "lgRMatrix")){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL(classes=extract_row_lgRMatrix(detection,i), parameters=parameters, reference.prob=ref$Q, pseudo=ref$pseudo)
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
        D_KL.randomized[r] <- get_D_KL(classes=classes,
                                       parameters=parameters, reference.prob=ref$Q, pseudo=ref$pseudo)
      }
      all.D_KL.randomized[i,] <- D_KL.randomized
    }# end for all T counts to select
  } else {

    for(i in 1:T.counts.to.select){
      setTxtProgressBar(pb, i) # progress bar

      T.count <- T.counts.selected[i]

      D_KL.randomized <- rep(NA,randomization.count)
      vector.to.randomize <- c(rep(TRUE, T.count), rep(FALSE, ncol(detection) - T.count))
      for(r in 1:randomization.count){
        # using default sampling
        D_KL.randomized[r] <- get_D_KL(classes=sample(x = vector.to.randomize),
                                       parameters=parameters, reference.prob=ref$Q, pseudo=ref$pseudo)
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

  message("### returning result...")
  # prepare the 'haystack' object to return
  res <- list(
    results = data.frame(
      D_KL = D_KL.observed,
      log.p.vals = p.vals,
      log.p.adj = p.adjs,
      T.counts = T.counts,
      row.names = row.names(detection)
    )
  )
  class(res) <- "haystack"
  res
}

#' Function to get the density of points with value TRUE in the (x,y) plot
#'
#' @param x x-axis coordinates of cells in a 2D representation (e.g. resulting from PCA or t-SNE)
#' @param y y-axis coordinates of cells in a 2D representation
#' @param detection A logical matrix or dgRMatrix showing which gens (rows) are detected in which cells (columns)
#' @param rows.subset Indices of the rows of 'detection' for which to get the densities. Default: all.
#' @param high.resolution Logical: should high resolution be used? Default is FALSE.
#'
#' @return A 3-dimensional array (dim 1: genes/rows of expression, dim 2 and 3: x and y grid points) with density data
get_density = function(x, y, detection, rows.subset=1:nrow(detection), high.resolution = FALSE){

  # set the parameters for getting the densities
  parameters <- get_parameters_haystack(x,y,high.resolution)

  densities <- array(data=NA, dim=c(length(rows.subset),parameters$grid.points))

  cl <- T # we are only looking at the T points here

  # detection could be a matrix class object now,
  # or a lgRMatrix object
  # if detection is a lgCMatrix object at this point, something is wrong
  if(is.matrix(detection)){
    for(i in 1:length(rows.subset)){
      r <- rows.subset[i]
      x.subset <- detection[r,]==cl
      y.subset <- detection[r,]==cl
      density <- kde2d_faster(dens.x=parameters$dens.x[,x.subset],
                              dens.y=parameters$dens.y[,y.subset])
      densities[i,,] <- density / sum(density)
    }
  } else if(inherits(detection, "lgRMatrix")){
    for(i in 1:length(rows.subset)){
      r <- rows.subset[i]
      x.subset <- extract_row_lgRMatrix(detection,r)==cl
      y.subset <- extract_row_lgRMatrix(detection,r)==cl
      density <- kde2d_faster(dens.x=parameters$dens.x[,x.subset],
                              dens.y=parameters$dens.y[,y.subset])
      densities[i,,] <- density / sum(density)
    }
  } else {
    stop("'detection' must be a matrix or lgRMatrix")
  }
  # set dimension names to genes, and grid points of x and y axes
  dimnames(densities) <- list(rownames(detection)[rows.subset],
                         seq(parameters$limits[1],parameters$limits[2],length.out = parameters$grid.points[1]), # x grid
                         seq(parameters$limits[3],parameters$limits[4],length.out = parameters$grid.points[2])  # y grid
                         )
  densities
}

#' Shows the results of the 'haystack' analysis in various ways, sorted by significance. Priority of params is genes > p.value.threshold > n.
#'
#' @param res.haystack A 'haystack' result variable
#' @param n If defined, the top "n" sigificant genes will be returned. Default: NA, which shows all results.
#' @param p.value.threshold If defined, genes passing this p-value threshold will be returned.
#' @param gene If defined, the results of this (these) gene(s) will be returned.
#'
#' @return A table with a sorted subset of the 'haystack' result according to input parameters.
#' @export
#'
#' @examples
#' # using the toy example of the singleCellHaystack package
#' # define a logical matrix with detection of each gene (rows) in each cell (columns)
#' dat.detection <- dat.expression > 1
#'
#' # running haystack in default mode
#' res <- haystack(dat.tsne, detection=dat.detection, method = "2D")
#'
#' # below are variations for showing the results in a table
#' # 1. list top 10 biased genes
#' show_result_haystack(res.haystack = res, n =10)
#' # 2. list genes with p value below a certain threshold
#' show_result_haystack(res.haystack = res, p.value.threshold=1e-10)
#' # 3. list a set of specified genes
#' set <- c("gene_497","gene_386", "gene_275")
#' show_result_haystack(res.haystack = res, gene = set)
show_result_haystack = function(res.haystack, n=NA, p.value.threshold=NA, gene=NA){

  # check input
  if(missing(res.haystack))
    stop("Parameter 'res.haystack' ('haystack' result) is missing")
  if(class(res.haystack)!="haystack")
    stop("'res.haystack' must be of class 'haystack'")
  if(is.null(res.haystack$results))
    stop("Results seem to be missing from 'haystack' result. Is 'res.haystack' a valid 'haystack' result?")
  if(!is.na(n) & !is.numeric(n))
    stop("The value of 'n' should be an integer")
  if(!is.na(n) & n > nrow(res.haystack$results))
    warning("Integer value of 'n' is larger than the number of rows in the 'haystack' results. Will return all results sorted.")
  if(!is.na(p.value.threshold) & (p.value.threshold<0 | p.value.threshold>1))
    stop("If 'p.value.threshold' is given as input, it should be between 0 and 1")

  # run through filters, one by one, if they are not NA
  # priority is genes > p.value.threshold > n
  result <- res.haystack$results
  if(any(!is.na(gene))){
    result <- result[is.element(rownames(result), gene), ]
  }
  if(!is.na(p.value.threshold)){
    result <- result[result[["log.p.vals"]] <= log10(p.value.threshold), ]
  }

  # at this point: 1) decide no. to return, and 2) sort by significance
  n.to.select <- ifelse(is.na(n), nrow(result), min(n, nrow(result)))
  o <- order(result$log.p.vals)
  result[o[1:n.to.select],]
}



