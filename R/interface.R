
#' Interface to various Haystack functions.
#'
#' @param object a matrix with expression values for genes in cells/spots, or a logical matrix indicating detection values, or a Seurat/SingleCellexperiment object containing expression values.
#' @param coordinates a matrix with coordinates of cells in PCA space (>=2 PCs) or other type of embedding (tSNE/UMAP), or 2D/3D spatial coordinates.
#' @param type highD, 2D, or spatial, depending on the type of input coordinates.
#' @param ... other parameters to pass on to the haystack functions.
#'
#' @param reduction name for the embeddings to be used when using Seurat or SingleCellExperiment objects.
#' @param assay name of the assay to be used when using Seurat or SingleCellExperiment objects.
#' @param slot name of the slot to be used when using Seurat objects.
#' @param dims dimensions to use.
#'
#' @return An object of class "haystack", including the results of the analysis, and the coordinates of the grid points used to estimate densities.
#' @export
#'
#' @examples
#' # I need to add some examples.
#' # A toy example will be added too.
haystack_interface = function(object, ...) {
  #stop("This function should not be used, or has mistakenly been marked for removal")

  # I think we can change the name of this function to just "haystack" or "singleCellHaystack"
  # we might have to do something about the s3.R function "haystack" in that case, to avoid confusion

  # might add several checks on the input here?
  # in that case, we can remove redundant, general checks from all the functions
  # and keep only specific checks inside the individual functions
  UseMethod("haystack_interface")
}

#' @rdname haystack_interface
#' @export
haystack_interface.matrix = function(object = NULL, coordinates = NULL, type = NULL, ...){
  #stop("This function should not be used, or has mistakenly been marked for removal")

  haystack_interface_raw(expression = object, coordinates = coordinates, type = type, ...)
}

#' @rdname haystack_interface
#' @export
haystack_interface.Seurat = function(object = NULL, reduction = "pca", assay = NULL, slot = "data", type = NULL, dims = NULL, ...) {
  #stop("This function should not be used, or has mistakenly been marked for removal")

  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package \"SeuratObject\" needed for this function to work. Please install it.", call. = FALSE)
  }

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(object)
  }
  expression <- SeuratObject::GetAssayData(object, slot = slot, assay = assay)

  coordinates <- SeuratObject::Embeddings(object, reduction)
  if (!is.null(dims)) {
    coordinates <- coordinates[, dims]
  }

  haystack_interface_raw(expression = expression, coordinates = coordinates, type = type, ...)
}


haystack_interface_raw <- function(expression = NULL, coordinates = NULL, type = NULL, ...) {
  #stop("This function should not be used, or has mistakenly been marked for removal")

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
    return(NULL)
  } else {
    # binary
    if (is.logical(expression)) {
      if (type == "highD") {
        # highD binary
        message("### calling haystack version highD binary...")
        res <- haystack_highD(x = coordinates, detection = expression, ...)
      }

      if (type == "2D") {
        # 2D binary
        message("### calling haystack version 2D binary...")
        res <- haystack_2D(x = coordinates[, 1], y = coordinates[, 2], detection = expression, ...)
      }

      if(type == "spatial") {
        # spatial binary
        # sorry: no such function!
        # and no plan to make it, really
        message("### there is no function to analyze spatial data using non-continuous (binary) expression values...")
      }
    }

    # continuous
    if (is.numeric(expression) || is(expression, "Matrix")) {
      if (type == "highD") {
        # highD continuous
        message("### calling haystack version highD continuous...")
        res <- haystack_continuous_highD(x = coordinates, expression = expression, ...)
      }
    }
  }

  res
}


