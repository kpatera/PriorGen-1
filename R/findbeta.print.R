#' The findbeta plot function
#'
#' A function that prints a summary any object of the class PriorGen
#'
#' @rdname print.PriorGen
#' @param x An object of type findbeta produced of one of the main PriorGen functions.
#' @examples
#' ## Example 1
#' ## Based on the available literature the mean value for the sensitivity of a test
#' ## is expected to be generally low and its variance not that low but not that much neither.
#' res_abs_1 <- findbeta_abstract(themean.cat = "Low", thevariance.cat = "Average")
#' print(res_abs_1)
#'
#' ## Example 2
#' ## Hierarchical prior
#' res_mult_1 <- findbetamupsi(
#'   themean = 0.10, percentile = 0.79,
#'   lower.v = TRUE, percentile.value = 0.26, psi.percentile = 0.95,
#'   percentile.median = 0.28, percentile95value = 0.3
#' )
#' print(res_mult_1)
#'
#' @export
print.PriorGen <- function(x) {
  findbeta.object <- x
  classchk <- as.character(class(findbeta.object))
  if (classchk != "PriorGen" && classchk != "PriorGen2") {
    stop("Provide an object of class PriorGen" & class(findbeta.object) != "PriorGen2")
  }
  if (length(findbeta.object) == 3) {
    paste0(
      "The desired Beta distribution that satisfies the specified conditions is: Beta(",
      round(findbeta.object$parameters[1], 2), ",", round(findbeta.object$parameters[2], 2), ").",
      " Verification: The percentile value ",
      round(qbeta(findbeta.object$input[2], findbeta.object$parameters[1], findbeta.object$parameters[2]), 2),
      " corresponds to the ", findbeta.object$input[2]*100, "th percentile"
    )

    # cat(paste("Descriptive statistics for this distribution can be found below:\n",findbeta.object$summary))
  } else {
    outgam <- as.vector(as.numeric(findbeta.object$param_gamma))
    paste0(
      "The desired Beta distribution that satisfies the specified conditions about the mean of the prevalence 'mu' is: Beta(",
      round(findbeta.object$param_beta[1], 2), round(findbeta.object$param_beta[2], 2),
      ") and The desired Gamma distribution that satisfies the specified conditions about the variability 'psi' of the prevalence is: Gamma(",
      as.numeric(round(outgam[1], 2)), as.numeric(round(outgam[2], 2)), ")"
    )
    # cat("Descriptive statistics for this distribution are:", findbeta.object$summary)
  }
  # cat(findbeta.object)
}
