

### write_haystack
#' Function to write haystack result data to file.
#'
#' @param res.haystack A 'haystack' result variable
#' @param file A file to write to
#'
#' @export
write_haystack = function (res.haystack, file){
  .Deprecated(msg = "This function has been deprecated and will be removed in the future.")

  # check input
  if(missing(res.haystack))
    stop("Parameter 'res.haystack' ('haystack' result) is missing")
  if(! inherits(res.haystack, "haystack"))
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
read_haystack = function (file){
  .Deprecated(msg = "This function has been deprecated and will be removed in the future.")

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

