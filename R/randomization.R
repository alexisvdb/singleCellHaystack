#' plot_rand_fit
#'
#' @param x haystack object.
#' @param type whether to plot mean or sd.
#'
#' @export
#'
plot_rand_fit <- function(x, type=c("mean", "sd")) {
  UseMethod("plot_rand_fit")
}

#' @rdname plot_rand_fit
#' @export
plot_rand_fit.haystack <- function(x, type=c("mean", "sd")) {
  type <- match.arg(type)

  d <- x$info$randomization[[type]]

  degree <- d$cv.selected["degree"]
  df <- d$cv.selected["df"]
  rmsd <- format(d$cv.selected["rmsd"], digits=2)

  spline_method <- x$info$randomization$method

  if (spline_method == "bs")
    subtitle <- paste0("spline=", spline_method, " degree=", degree, " df=", df, " rmsd=", rmsd)

  if (spline_method == "ns")
    subtitle <- paste0("spline=", spline_method, " df=", df, " rmsd=", rmsd)

  method <- strsplit(x$info$method, "_")[[1]]

  if (method[1] == "continuous")
    xlab <- "logCV"

  if (method[1] == "binary")
    xlab <- "# positive cells"

  ylab <- paste0(type, " of logKLD")

  ggplot(d$observed, aes(.data[["x"]], .data[["y"]])) +
    geom_point() +
    geom_line(color="red", data=d$fitted) +
    labs(
      x = xlab,
      y = ylab,
      title=paste0(x$info$method, ": fitted values for ", type),
      subtitle=subtitle)
}
