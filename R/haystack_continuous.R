
#' The main Haystack function, for higher-dimensional spaces and continuous expression levels.
#'
#' @param x Coordinates of cells in a 2D or higher-dimensional space. Rows represent cells, columns the dimensions of the space.
#' @param expression a matrix with expression data of genes (rows) in cells (columns)
#' @param grid.points An integer specifying the number of centers (grid points) to be used for estimating the density distributions of cells. Default is set to 100.
#' @param weights.advanced.Q (Default: NULL) Optional weights of cells for calculating a weighted distribution of expression.
#' @param dir.randomization If NULL, no output is made about the random sampling step. If not NULL, files related to the randomizations are printed to this directory.
#' @param scale Logical (default=TRUE) indicating whether input coordinates in x should be scaled to mean 0 and standard deviation 1.
#' @param grid.method The method to decide grid points for estimating the density in the high-dimensional space. Should be "centroid" (default) or "seeding".
#' @param randomization.count Number of randomizations to use. Default: 100
#' @param n.genes.to.randomize Number of genes to use in randomizations. Default: 100
#' @param selection.method.genes.to.randomize Method used to select genes for randomization.
#' @param grid.coord matrix of grid coordinates.
#' @param spline.method Method to use for fitting splines "ns" (default): natural splines, "bs": B-splines.
#'
#' @return An object of class "haystack", including the results of the analysis, and the coordinates of the grid points used to estimate densities.
#' @export
#'
#' @examples
#' # using the toy example of the singleCellHaystack package
#'
#' # running haystack
#' res <- haystack(dat.tsne, dat.expression)
#' # list top 10 biased genes
#' show_result_haystack(res, n=10)
haystack_continuous_highD = function(x, expression, grid.points = 100, weights.advanced.Q=NULL,
                                     dir.randomization = NULL, scale=TRUE, grid.method="centroid",
                                     randomization.count = 100, n.genes.to.randomize = 100,
                                     selection.method.genes.to.randomize = "heavytails", grid.coord=NULL,
                                     spline.method="ns"){

  renderPB <- isFALSE(getOption("rstudio.notebook.executing"))

  message("### calling haystack_continuous_highD()...")

  if (!is.null(grid.coord)) {
    grid.points <- nrow(grid.coord)

    if (ncol(x) != ncol(grid.coord))
      stop("Coordinates and grid points have different number of columns.")

    if (!is.null(colnames(x)) || !is.null(colnames(grid.coord)))
      if (! identical(colnames(x), colnames(grid.coord)))
        stop("Coordinates and grid points have different column names.")
  }

  # Check for sparseMatrixStats package.
  if(requireNamespace("sparseMatrixStats", quietly = TRUE)){
    message("### Using package sparseMatrixStats to speed up statistics in sparse matrices.")
    useSMS <- TRUE
  } else {
    message("### Package sparseMatrixStats not found. Install for speed improvements.")
    useSMS <- FALSE
  }

  message("### Calculating row-wise mean and SD... ")
  if (useSMS) {
    #expr.mean <- sparseMatrixStats::rowMeans2(expression)
    #names(expr.mean) <- rownames(expression)
    expr.mean <- Matrix::rowMeans(expression)
  }
  else {
    #expr.mean <- apply(expression,1,mean)
    expr.mean <- Matrix::rowMeans(expression) # Maybe this is best as default?
  }

  if (min(expr.mean) < 0)
    stop("Some features have an average signal < 0. Expect average signal >= 0.")

  if (useSMS)
    expr.sd <- sparseMatrixStats::rowSds(expression)
  else
    expr.sd <- apply(expression,1,sd)

  sel.bad <- expr.sd == 0
  expression <- expression[!sel.bad, ]
  expr.mean <- expr.mean[!sel.bad]
  expr.sd <- expr.sd[!sel.bad]
  message("### Filtered ", sum(sel.bad)," genes with zero variance...")

  expr.mean <- expr.mean + 1e-300
  expr.sd <- expr.sd + 1e-300

  message("### Using ",randomization.count," randomizations...")
  message("### Using ",n.genes.to.randomize," genes to randomize...")

  # check input
  if(!is.numeric(x) && !all(apply(x, 2, is.numeric)))
    stop("'x' must be a numeric matrix or data.frame")
  if(!is.matrix(x) && !is.data.frame(x))
    stop("'x' must be a numeric matrix or data.frame")
  #if(ncol(x) < 2)
  #  stop("'x' must have at least 2 columns")
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
  if(spline.method!="ns" & spline.method!="bs")
    stop("'spline.method' should be either 'ns' or 'bs'")

  # if expression is a dgCMatrix, convert it to a dgRMatrix
  if(inherits(expression, "dgCMatrix")){
    message("### converting expression data from dgCMatrix to dgRMatrix")
    expression <- as(expression, "RsparseMatrix")
  }

  count.cells <- ncol(expression)
  count.genes <- nrow(expression)
  if(n.genes.to.randomize > count.genes){
    warning("Number of genes to randomize (",n.genes.to.randomize,") is higher than the number of genes in the data.")
    warning("Setting number of genes to randomize to ",count.genes)
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
  if (is.null(grid.coord)) {
    message("### deciding grid points...")
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
  bandwidth <- median(apply(dist.to.grid,1,min))
  dist.to.grid.norm <- dist.to.grid / bandwidth
  density.contributions <-
    exp(-dist.to.grid.norm * dist.to.grid.norm / 2)

  if(is.null(weights.advanced.Q)){
    Q <- colSums(density.contributions)
  } else {
    Q <- colSums(density.contributions*weights.advanced.Q)
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
  if (renderPB)
    pb <- txtProgressBar(min = 0, max = count.genes, style = 3, file = stderr()) # progress bar
  if(is.matrix(expression)){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL_continuous_highD(weights=expression[i,], density.contributions = density.contributions, reference.prob = Q, pseudo = pseudo)
      if (renderPB) setTxtProgressBar(pb, i) # progress bar
    }
  } else if(inherits(expression, "dgRMatrix")){
    for(i in 1:count.genes){
      D_KL.observed[i] <- get_D_KL_continuous_highD_SPARSE(weights_list=extract_row_dgRMatrix_as_sparse(expression,i), density.contributions = density.contributions, reference.prob = Q, pseudo = pseudo)
      if (renderPB) setTxtProgressBar(pb, i) # progress bar
    }
  }
  if (renderPB) close(pb) # progress bar

  ### use randomization to estimate expected values
  ### and p values for D_KL
  # for the mean D_KL:  the mean(log(D_KL)) can be modeled in function of log(coeffVar)
  # for the SD of D_KL: the sd(log(D_KL))   can be modeled in function of coeffVar
  # we need the Coefficient of Variation (coeffVar) ---> from the input data (expression)
  # we need the D_KL values of (a small set of) randomized genes ---> do randomizations

  coeffVar  <- expr.sd/expr.mean # coefficient of variation

  message("### performing randomizations...")

  # use all genes if the total number of genes is less than n.genes.to.randomize
  o <- order(coeffVar)
  # Option "uniform"   : evenly spread ("uniform")
  # Option "heavytails": top 10 and bottom 10, otherwise evenly spread
  if(selection.method.genes.to.randomize=="uniform"){
    genes.to.randomize <- o[floor(seq(1,count.genes, length.out = n.genes.to.randomize))]
  } else if(selection.method.genes.to.randomize=="heavytails"){
    genes.to.randomize <- o[c(1:9,
                            as.integer(seq(10,count.genes-9,length.out = n.genes.to.randomize-18)),
                            length(o)-(8:0))]
  }


  # - do x randomizations and get their D_KL values
  # - get mean and SD of D_KL values,
  # - do fitting in function of coeffVar
  # - use those results to estimate p values of all genes

  # to store randomization results
  all.D_KL.randomized <- matrix(NA,nrow=n.genes.to.randomize, ncol=randomization.count)

  if (renderPB) pb <- txtProgressBar(min = 0, max = n.genes.to.randomize, style = 3, file = stderr()) # progress bar


  if(is.matrix(expression)){
    for(i in 1:n.genes.to.randomize){
      if (renderPB) setTxtProgressBar(pb, i) # progress bar

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
  } else if(inherits(expression, "dgRMatrix")){
    for(i in 1:n.genes.to.randomize){
      if (renderPB) setTxtProgressBar(pb, i) # progress bar

      D_KL.randomized <- rep(NA,randomization.count)
      weights.list.to.randomize <- extract_row_dgRMatrix_as_sparse(expression,genes.to.randomize[i])
      non.zero.n <- length(weights.list.to.randomize$ind)
      for(r in 1:randomization.count){
        # using default sampling
        weights.list.to.randomize$ind <- sample(x = count.cells, size = non.zero.n)
        D_KL.randomized[r] <- get_D_KL_continuous_highD_SPARSE(
          weights_list=weights.list.to.randomize,
          density.contributions = density.contributions, reference.prob = Q, pseudo = pseudo
        )
      }
      all.D_KL.randomized[i,] <- D_KL.randomized
    }# end for all T counts to select
  }
  if (renderPB) close(pb) # progress bar


  message("### estimating p-values...")
  p.vals <- get_log_p_D_KL_continuous(D_KL.observed = D_KL.observed,
                                      D_KL.randomized = all.D_KL.randomized,
                                      all.coeffVar = coeffVar,
                                      train.coeffVar = coeffVar[genes.to.randomize],
                                      output.dir = dir.randomization,
                                      spline.method = spline.method)
  rand_info <- p.vals$info
  rand_info$method <- spline.method
  rand_info$genes_to_randomize <- genes.to.randomize
  p.vals <- p.vals$fitted

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
    info = list(
      method="continuous_highD",
      randomization = rand_info,
      grid.coordinates = grid.coord,
      coord_mean = x.scale.center,
      coord_std = x.scale.scale,
      densities = density.contributions,
      cv = coeffVar
    )
     #,
    #all.D_KL.randomized = all.D_KL.randomized
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
  P <- colSums(density.contributions*weights)
  P <- P + pseudo
  P <- P / sum(P)
  D_KL <- sum(P * log(P/Q))

  D_KL
}


#' Estimates the significance of the observed Kullback-Leibler divergence by comparing to randomizations for the continuous version of haystack.
#'
#' @param D_KL.observed A vector of observed Kullback-Leibler divergences.
#' @param D_KL.randomized A matrix of Kullback-Leibler divergences of randomized datasets.
#' @param all.coeffVar Coefficients of variation of all genes. Used for fitting the Kullback-Leibler divergences.
#' @param train.coeffVar Coefficients of variation of genes that will be used for fitting the Kullback-Leibler divergences.
#' @param output.dir Optional parameter. Default is NULL. If not NULL, some files will be written to this directory.
#' @param spline.method Method to use for fitting splines "ns" (default): natural splines, "bs": B-splines.
#'
#' @return A vector of log10 p values, not corrected for multiple testing using the Bonferroni correction.
get_log_p_D_KL_continuous = function(D_KL.observed, D_KL.randomized, all.coeffVar, train.coeffVar, output.dir = NULL, spline.method="ns"){

  # prepare data for fitting
  log.D_KL <- log(D_KL.randomized)
  D_KL.log.mean <- apply(log.D_KL,1,mean)
  D_KL.log.sd   <- apply(log.D_KL,1,sd)


  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  # first, for the mean D_KL
  message("### picking model for mean D_KL...")
  plot.file <- NULL
  if(!is.null(output.dir))
    plot.file <- paste0(output.dir,"/fit_logCoeffVar_vs_meanLogD_KL.pdf")
  if(spline.method=="ns"){
    model <- get_model_cv_ns(x = log(train.coeffVar), y = D_KL.log.mean, plot.file = plot.file)
  } else if(spline.method=="bs"){
    model <- get_model_cv_bs(x = log(train.coeffVar), y = D_KL.log.mean, plot.file = plot.file)
  }

  info <- list()
  info$features <- rownames(D_KL.randomized)
  info$mean <- model$info
  model <- model$model
  #summary(model)

  # fitted values for all cases
  fitted.D_KL_log.mean <- predict(model, data.frame(x=log(all.coeffVar)), type="response")


  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  # second, for the SD of the D_KL
  message("### picking model for stdev D_KL...")
  plot.file <- NULL
  if(!is.null(output.dir))
    plot.file <- paste0(output.dir,"/fit_logCoeffVar_vs_SdLogD_KL.pdf")
  if(spline.method=="ns"){
    model <- get_model_cv_ns(x = log(train.coeffVar), y = D_KL.log.sd, plot.file = plot.file)
  } else if(spline.method=="bs"){
    model <- get_model_cv_bs(x = log(train.coeffVar), y = D_KL.log.sd, plot.file = plot.file)
  }
  info$sd <- model$info
  model <- model$model
  #summary(model)

  # fitted values for all cases
  fitted.D_KL_log.sd <- predict(model, data.frame(x=log(all.coeffVar)), type="response")


  # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  # estimate p values
  fitted.log.p.vals <- pnorm(log(D_KL.observed), mean = fitted.D_KL_log.mean, sd = fitted.D_KL_log.sd, lower.tail = FALSE, log.p = T)/log(10)

  list(fitted=fitted.log.p.vals, info=info)
}


get_model_cv_ns = function(x, y, plot.file = NULL){
  message("### using natural splines")
  # prepare cross validation
  n_cv <- 10 # 10 fold
  sets <- sample(rep(1:n_cv, length.out = length(x))) # I am using sample() to really randomize them

  # set boundary knots to be slightly outside range of predictor values
  range.x <- range(x)
  r.x <- range.x[2]-range.x[1]
  boundary.knots <- c(range.x[1] - 0.01*r.x, range.x[2] + 0.01*r.x)

  dfs <- 1:10

  best_rmsd   <- Inf
  best_df     <- NA
  for(df in dfs){

    rmsds <- rep(NA, n_cv)
    for(s in 1:n_cv){
      x.train <- x[sets!=s]
      y.train <- y[sets!=s]
      x.test  <- x[sets==s]
      y.test  <- y[sets==s]

      model <- lm(y.train ~ ns(x.train, df=df,Boundary.knots=boundary.knots))

      y.pred <- predict(model, data.frame(x.train=x.test), type="response")
      rmsds[s] <- sqrt(mean((y.pred - y.test)^2))
    }
    rmsd <- mean(rmsds)

    # message(degree," ",df," ",rmsd)

    # update if better solution
    if(rmsd < best_rmsd){
      best_rmsd <- rmsd
      best_df <- df
    }
  }# for dfs


  message("### best RMSD  : ",round(best_rmsd,3))
  message("### best df    : ",best_df)

  model <- lm(y ~ ns(x, df=best_df,Boundary.knots=boundary.knots))

  x.seq <- seq(range.x[1],range.x[2],length.out=1000)
  fitted.y <- predict(model, data.frame(x=x.seq), type="response")

  info <- list(
    cv.parameters=list(dfs=dfs),
    cv.selected=list(df=best_df, rmsd=best_rmsd),
    observed=data.frame(feature=names(x), x=x, y=y),
    fitted=data.frame(feature=names(x), x=x.seq, y=fitted.y),
    spline.method = "ns"
  )

  if(!is.null(plot.file)){
    range.y <- range(range(y), range(fitted.y))
    pdf(plot.file)
    plot(x, y, ylim = range.y)
    points(x.seq,fitted.y,col="red", type="l")
    dev.off()
  }

  list(model=model, info=info)
}

get_model_cv_bs = function(x, y, plot.file = NULL){
  message("### using B-splines")
  # prepare cross validation
  n_cv <- 10 # 10 fold
  sets <- sample(rep(1:n_cv, length.out = length(x))) # I am using sample() to really randomize them

  # set boundary knots to be slightly outside range of predictor values
  range.x <- range(x)
  r.x <- range.x[2]-range.x[1]
  boundary.knots <- c(range.x[1] - 0.01*r.x, range.x[2] + 0.01*r.x)

  degrees <- 1:5
  dfs <- 1:10

  best_rmsd   <- Inf
  best_degree <- NA
  best_df     <- NA
  for(degree in degrees){
    for(df in dfs){

      # check if df is suitable
      # see: https://github.com/SurajGupta/r-source/blob/master/src/library/splines/R/splines.R
      ord <- 1L + (degree <- as.integer(degree))
      intercept <- FALSE
      nIknots <- df - ord + (1L - intercept)
      if(nIknots < 0L){
        # message("degree: ",degree, " - df: ",df," - df too small; skip")
        next
      }

      rmsds <- rep(NA, n_cv)
      for(s in 1:n_cv){
        x.train <- x[sets!=s]
        y.train <- y[sets!=s]
        x.test  <- x[sets==s]
        y.test  <- y[sets==s]

        model <- lm(y.train ~ bs(x.train, df=df,Boundary.knots=boundary.knots, degree=degree))

        y.pred <- predict(model, data.frame(x.train=x.test), type="response")
        rmsds[s] <- sqrt(mean((y.pred - y.test)^2))
      }
      rmsd <- mean(rmsds)

      # message(degree," ",df," ",rmsd)

      # update if better solution
      if(rmsd < best_rmsd){
        best_rmsd <- rmsd
        best_degree <- degree
        best_df <- df
      }
    }# for dfs
  }# for degrees

  message("### best RMSD  : ",round(best_rmsd,3))
  message("### best degree: ",best_degree)
  message("### best df    : ",best_df)

  model <- lm(y ~ bs(x, df=best_df,Boundary.knots=boundary.knots, degree=best_degree))

  x.seq <- seq(range.x[1],range.x[2],length.out=1000)
  fitted.y <- predict(model, data.frame(x=x.seq), type="response")

  info <- list(
    cv.parameters=list(degrees=degrees, dfs=dfs),
    cv.selected=list(degree=best_degree, df=best_df, rmsd=best_rmsd),
    observed=data.frame(feature=names(x), x=x, y=y),
    fitted=data.frame(feature=1:length(x.seq), x=x.seq, y=fitted.y),
    spline.method = "bs"
  )

  if(!is.null(plot.file)){
    range.y <- range(range(y), range(fitted.y))
    pdf(plot.file)
    plot(x, y, ylim = range.y)
    points(x.seq,fitted.y,col="red", type="l")
    dev.off()
  }

  list(model=model, info=info)
}



get_D_KL_continuous_highD_SPARSE = function(weights_list, density.contributions, reference.prob, pseudo = 0){

  val = weights_list$val
  ind = weights_list$ind

  # the reference distribution Q of cells in the space
  Q <- reference.prob

  # calculating the Kullback-Leibler divergence of the distribution
  # of expression vs reference distribution Q
  P <- colSums(density.contributions[ind,,drop=FALSE]*val)
  P <- P + pseudo
  P <- P / sum(P)
  D_KL <- sum(P * log(P/Q))

  D_KL
}

prepare_clustering <- function(x, expression, grid.coordinates, scale=TRUE) {
  ngenes <- nrow(expression)
  # scale data if needed
  if(scale){
    x <- scale(x)
    # save the mean and stdev of the scaling
    x.scale.center <- attr(x = x, which = "scaled:center")
    x.scale.scale <- attr(x = x, which = "scaled:scale")
    grid.coordinates <- (grid.coordinates - rep(x.scale.center,each=nrow(grid.coordinates))) / rep(x.scale.scale,each=nrow(grid.coordinates))
  }

  # if expression is a dgCMatrix, convert it to a dgRMatrix
  if(inherits(expression, "dgCMatrix")){
    message("### converting expression data from dgCMatrix to dgRMatrix")
    expression <- as(expression, "RsparseMatrix")
  }

  # get densities
  #expression.rownames <- rownames(expression)
  #row.index.subset <- which(is.element(expression.rownames, genes))

  dist.to.grid <- get_dist_two_sets(x,grid.coordinates)

  # process the distances to a suitable density contribution
  # first, set bandwidth
  # bandwidth <- sqrt(sum((apply(x, 2, default_bandwidth.nrd)) ^ 2))
  bandwidth <- median(apply(dist.to.grid,1,min))
  dist.to.grid.norm <- dist.to.grid / bandwidth
  density.contributions <-
    exp(-dist.to.grid.norm * dist.to.grid.norm / 2)

  #densities <- matrix(NA, nrow=length(row.index.subset), ncol=ncol(density.contributions))
  densities <- matrix(NA, nrow=ngenes, ncol=ncol(density.contributions))
  rownames(densities) <- rownames(expression)
  #row.names(densities) <- expression.rownames[row.index.subset]

  message("### collecting density data...")
  pb <- txtProgressBar(min = 0, max = ngenes, style = 3, file = stderr()) # progress bar
  if(is.matrix(expression)){
    #for(g in 1:length(row.index.subset)){
    for(g in seq_len(ngenes)){
      #gene_index <- row.index.subset[g]
      #densities[g,] <- colSums(density.contributions*expression[gene_index,])
      densities[g,] <- colSums(density.contributions*expression[g,])
      setTxtProgressBar(pb, g) # progress bar
    }
  } else if( inherits(expression, "dgRMatrix") ){
    #for(g in 1:length(row.index.subset)){
    for(g in seq_len(ngenes)){
      #gene_index <- row.index.subset[g]
      #densities[g,] <- colSums(density.contributions*extract_row_dgRMatrix(expression,gene_index))
      densities[g,] <- colSums(density.contributions*extract_row_dgRMatrix(expression,g))
      setTxtProgressBar(pb, g) # progress bar
    }
  } else {
    stop("'expression' must be a matrix or dgRMatrix")
  }
  close(pb) # progress bar

  # rescale to sum to 1. This is to avoid R thinking sd=0 in the case where an entire row has very low values
  densities <- densities / rowSums(densities)
  densities
}

hclust_haystack_continuous = function(x, expression, grid.coordinates, hclust.method="ward.D", cor.method="spearman", scale = TRUE, ...){

  # many checks on input
  # see haystack_continuous_highD and also
  # see the hclust_haystack_highD function

  densities <- prepare_clustering(x, expression, grid.coordinates, scale=scale)

  #heatmap(dist.to.grid.norm, Rowv=NA, Colv=NA, scale="none")
  #heatmap(densities, Rowv=NA, Colv=NA, scale="none")

  dist <- as.dist(1 - cor(t(densities),method=cor.method)) # dist(densities)
  hc <- hclust(dist, method=hclust.method, ...)
  hc
}

kmeans_haystack_continuous = function(x, expression, grid.coordinates, k, scale = TRUE, ...){

  # many checks on input
  # see haystack_continuous_highD and also
  # see the kmeans_haystack_highD function

  densities <- prepare_clustering(x, expression, grid.coordinates, scale=scale)

  #heatmap(dist.to.grid.norm, Rowv=NA, Colv=NA, scale="none")
  #heatmap(densities, Rowv=NA, Colv=NA, scale="none")

  km <- kmeans(x=densities, centers=k, ...)
  km
}
