#' The findbeta lines function
#'
#' A function that plots another PriorGen object on top of a class PriorGen plot.
#'
#' @rdname lines.PriorGen
#' @param x An object of type PriorGen produced of 
#' one of the other findbeta functions.
#' @param ... More basic plot arguments
#' @examples
#' ## Example 1
#' ## Based on the available literature the mean value
#' ## for the sensitivity of a test is expected to be
#' ## generally low and its variance not that low but not that much neither.
#'
#' res_abs_1 <- findbeta_abstract(
#'   themean.cat = "Low",
#'   thevariance.cat = "Average"
#' )
#'
#' plot(res_abs_1,
#'   main = "Plot of the findbeta_abstract function",
#'   lwd = 3, ylim = c(0, 7), type = "l"
#' )
#'
#' res_abs_2 <- findbeta_abstract(
#'   themean.cat = "High",
#'   thevariance.cat = "Average"
#' )
#'
#' lines(res_abs_2, lwd = 3, col = "red")
#'
#' ## Example 2
#' ## Hierarchical prior
#'
#' res_mult_1 <- findbetamupsi(
#'   themean = 0.10, percentile = 0.79,
#'   lower.v = TRUE, percentile.value = 0.26, psi.percentile = 0.95,
#'   percentile.median = 0.28, percentile95value = 0.3
#' )
#'
#' plot(res_mult_1,
#'   main = "Plot of the findbetamupsi function",
#'   lwd = 3, ylim = c(0, 7)
#' )
#'
#' res_mult_2 <- findbetamupsi(
#'   themean = 0.12, percentile = 0.79,
#'   lower.v = TRUE, percentile.value = 0.26, psi.percentile = 0.95,
#'   percentile.median = 0.28, percentile95value = 0.3
#' )
#'
#' lines(res_mult_2, lwd = 3, col = "red")
#' @export


lines.PriorGen <- function(x, ...) {
  findbeta.object <- x
  classchk <- as.character(class(findbeta.object))
  if (classchk != "PriorGen" && classchk != "PriorGen2") {
    stop("Provide an object of class PriorGen")
  }

  if (length(findbeta.object) == 3) {
    a_plot <- findbeta.object$parameters[1]
    b_plot <- findbeta.object$parameters[2]
    x <- seq(0, 1, 0.001)
    y <- dbeta(seq(0, 1, 0.001), a_plot, b_plot)
    lines(x, y, ...)
  } else {
    a_plot <- mean(findbeta.object$param_upper$at)
    b_plot <- mean(findbeta.object$param_upper$bt)
    x <- seq(0, 1, 0.001)
    y <- dbeta(seq(0, 1, 0.001), a_plot, b_plot)
    lines(x, y, ...)
  }
}
