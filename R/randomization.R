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

  d <- x$randomization[[type]]

  ggplot(d$observed, aes(x, y)) +
    geom_point() +
    geom_line(color="red", data=d$fitted)
}
