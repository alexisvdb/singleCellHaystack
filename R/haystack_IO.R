########################################
########################################
### write_haystack
#' Function to write haystack result data to file.
#'
#' @param x A 'haystack' result variable
#' @param file A file to write to
#'
#' @export
#'
#' @examples
#' warn("I will add this later")
write_haystack = function (x, file){

  # check input
  if(missing(x))
    stop("Parameter 'x' ('haystack' result) is missing")
  if(class(x)!="haystack")
    stop("'x' must be of class 'haystack'")
  if(is.null(x$results))
    stop("Results seem to be missing from 'haystack' result. Is 'x' a valid 'haystack' result?")
  if(missing(file))
    stop("Parameter 'file' is missing")

  write.csv(x = x$results, file = file)
}


########################################
########################################
### read_haystack
#' Function to read haystack results from file.
#'
#' @param file A file containing 'haystack' results to read
#'
#' @return An object of class "haystack"
#' @export
#'
#' @examples
#' warn("I will add this later")
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

