
#' Interface to various Haystack functions.
#'
#' @param coordinates typically >2D coordinates of cells in the PC space (2D is also fine), or 2D/3D spatial coordinates.
#' @param expression matrix or data.frame with continuous expression values, or binary (TRUE/FALSE) detection values of genes in cells/spots.
#' @param type highD, 2D, or spatial, depending on the type of input coordinates.
#' @param ... other parameters to pass on to the haystack functions.
#'
#' @return An object of class "haystack", including the results of the analysis, and the coordinates of the grid points used to estimate densities.
#' @export
#'
#' @examples
#' # I need to add some examples.
#' # A toy example will be added too.
haystack_interface = function(coordinates, expression, type, ...) {
  UseMethod("haystack_interface")
}

#' @rdname haystack_interface
#' @export
haystack_interface.matrix = function(coordinates = NULL, expression = NULL, type = NULL, ...){
  # I think we can change the name of this function to just "haystack" or "singleCellHaystack"
  # we might have to do something about the s3.R function "haystack" in that case, to avoid confusion

  # might add several checks on the input here?
  # in that case, we can remove redundant, general checks from all the functions
  # and keep only specific checks inside the individual functions

  if (is.null(coordinates) || is.null(expression) || is.null(type)) {
    message("### usage:\n\n",
            "    haystack_interface(coordinates, [expression or detection data], type)\n\n",
            "    where:\n",
            "    - coordinates: typically >2D coordinates of cells in the PC space (2D is also fine), \n",
            "                   or 2D/3D spatial coordinates\n",
            "    - expression : matrix or data.frame with continuous expression values of genes in cells/spots\n",
            "    - detection  : matrix or data.frame with binary (TRUE or FALSE) detection values \n",
            "                   of genes in cells/spots\n",
            "    - type       : highD, 2D, or spatial, depending on the type of input coordinates\n\n")
  } else {
    # binary
    if (mode(expression) == "logical") {
      if (type == "highD") {
        # highD binary
        message("### calling haystack version highD binary...")
        haystack_highD(x = coordinates, detection = expression)
      }

      if (type == "2D") {
        # 2D binary
        message("### calling haystack version 2D binary...")
        haystack_2D(x = coordinates[, 1], y = coordinates[, 2], detection = expression)
      }

      if(type == "spatial") {
        # spatial binary
        # sorry: no such function!
        # and no plan to make it, really
        message("### there is no function to analyze spatial data using non-continuous (binary) expression values...")
      }
    }

    # continuous
    if (mode(expression) == "numeric") {
      if (type == "highD") {
        # highD continuous
        message("### calling haystack version highD continuous...")
        haystack_continuous_highD(x = coordinates, expression = expression)
      }

      if (type == "2D") {
        # 2D continuous
        message("### calling haystack version 2D continuous...")
        haystack_continuous_2D(x = coordinates[, 1], y = coordinates[, 2], expression = expression)
      }

      if(type == "spatial") {
        # spatial continuous
        message("### calling haystack version spatial continuous...")
        haystack_continuous_spatial(coordinates = coordinates, expression = expression)
      }
    }
  }
}

