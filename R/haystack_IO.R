

### write_haystack
#' Function to write haystack result data to file.
#'
#' @param res.haystack A 'haystack' result variable
#' @param file A file to write to
#'
#' @export
#'
#' @examples
#' # using the toy example of the singleCellHaystack package
#' # define a logical matrix with detection of each gene (rows) in each cell (columns)
#' dat.detection <- dat.expression > 1
#'
#' # running haystack in default mode
#' res <- haystack(x=dat.tsne$tSNE1, y=dat.tsne$tSNE2, detection=dat.detection)
#'
#' # write result to file outfile.csv
#' write_haystack(res.haystack = res, file = "outfile.csv")
#'
#' # read in result from file
#' res.copy <- read_haystack(file = "outfile.csv")
write_haystack = function (res.haystack, file){

  # check input
  if(missing(res.haystack))
    stop("Parameter 'res.haystack' ('haystack' result) is missing")
  if(class(res.haystack)!="haystack")
    stop("'res.haystack' must be of class 'haystack'")
  if(is.null(res.haystack$results))
    stop("Results seem to be missing from 'haystack' result. Is 'res.haystack' a valid 'haystack' result?")
  if(missing(file))
    stop("Parameter 'file' is missing")

  write.csv(x = res.haystack$results, file = file)
}




### read_haystack
#' Function to read haystack results from file.
#'
#' @param file A file containing 'haystack' results to read
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
#' res <- haystack(x=dat.tsne$tSNE1, y=dat.tsne$tSNE2, detection=dat.detection)
#'
#' # write result to file outfile.csv
#' write_haystack(res.haystack = res, file = "outfile.csv")
#'
#' # read in result from file
#' res.copy <- read_haystack(file = "outfile.csv")
read_haystack = function (file){

  # check input
  if(missing(file))
    stop("Parameter 'file' is missing")

  x <- read.csv(file = file, row.names=1)

  # prepare the 'haystack' object to return
  res <- list(
    results = x
  )
  class(res) <- "haystack"
  res
}

