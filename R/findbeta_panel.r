#' The findbeta (panel) function
#'
#'
#' A function to estimate the parameters alpha and beta of a Beta distribution based on the existing prior beliefs (data and/or expert opinion). Information is provided about the mean (or the median or the mode) and a corresponding scale metric, either the variance or the range of the parameter.
#'
#' @usage findbeta_panel(themean.vec=NULL, themedian.vec=NULL,
#'  themode.vec=NULL, seed=280385, nsims=10000)
#'
#' @param themean.vec specify your prior belief about the mean. It takes a vector of means, with values between 0 and 1. Not to be specified if a vector has been given for the median or the mode.
#' @param themedian.vec specify your prior belief about the median. It takes a vector of medians, with values between 0 and 1. Not to be specified if a vector has been given for the mean or the mode.
#' @param themode.vec specify your prior belief about the mode. It takes a vector of modes, with values between 0 and 1. Not to be specified if a vector has been given for the mean or the median.
#' @param seed A fixed seed for replication purposes.
#' @param nsims Number of simulations for the creation of various summary metrics of the elicited prior.
#'
#' @examples
#' ## Example 1
#' ## Based on the available literature the median/mean/mode value for the specificity of a
#' ## test is expected to be equal to 0.1, 0.2, 0.4, 0.04, 0.01, 0.5 based on opinions of 6 experts.
#'
#' resmed <- findbeta_panel(themedian.vec = c(0.1, 0.2, 0.4, 0.04, 0.01, 0.5))
#' resmed
#'
#' resmea <- findbeta_panel(themean.vec = c(0.1, 0.2, 0.4, 0.04, 0.01, 0.5))
#' resmea
#'
#' resmod <- findbeta_panel(themode.vec = c(0.1, 0.2, 0.4, 0.04, 0.01, 0.5))
#' resmod
#'
#' plot(resmed, lty = 1)
#' lines(resmea, lty = 2)
#' lines(resmod, lty = 3)
#' @export
#' @return parameters: The beta distribution parameters Beta(a,b)
#' @return summary: A basic summary of the elicited prior
#' @return input: The initial input value that produced the above prior.
#'
#' @references
#' Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2005): Estimation of diagnostic test sensitivity and specificity through Bayesian modeling. Preventive veterinary medicine, \bold{68}, 145--163.

findbeta_panel <- function(themean.vec = NULL, themedian.vec = NULL, themode.vec = NULL,
                           seed = 280385, nsims = 10000) {
  stopifnot((is.null(themean.vec) + is.null(themedian.vec) + is.null(themode.vec)) == 2)
  themedian <- themean <- themode <- NULL
  # Delta method for transposing log(p) exp(log(p))
  if (!is.null(themean.vec)) {
    thevariance <- var(themean.vec)
    themean <- mean(themean.vec)
    therange <- c(0, 1)
  } else if (!is.null(themedian.vec)) {
    themedian <- mean(themedian.vec)
    therange <- c(min(themedian.vec), max(themedian.vec))
    thevariance <- NULL
  } else if (!is.null(themode.vec)) {
    themode <- mean(themode.vec)
    therange <- c(min(themode.vec), max(themode.vec))
    thevariance <- NULL
  }

  stopifnot(((thevariance <= 0.5) + (therange[2] <= 1) + (therange[1] >= 0)) == 3)
  out <- NULL
  alpha <- 0.99995
  percentile <- 0.9999

  if (!is.null(themedian)) {
    name <- "themedian"
    value <- themedian
  }
  if (!is.null(themode)) {
    name <- "themode"
    value <- themode
  }
  alpha <- 0.99995
  if (!is.null(themean)) {
    if ((themean + qnorm(alpha) * thevariance) > 1) {
      ("The set variance resulted in percentile.value greater than 1. The percentile.value was set to 1")
    }
    if ((themean + qnorm(alpha) * thevariance) > 1) {
      thevariance <- (1 - themean) / qnorm(alpha)
    }
  }
  if (!is.null(themean)) {
    name <- "themean"
    value <- themean
    scalevalue <- thevariance
    percentile.value <- themean + qnorm(alpha) * thevariance
    percentile.value <- percentile.value - 0.001
  }
  if (!is.null(themedian)) {
    name <- "themedian"
    value <- themedian
    scalevalue <- therange[2] - therange[1]
    percentile.value <- therange[2]
    percentile.value <- percentile.value - 0.001
  }
  if (!is.null(themode)) {
    name <- "themode"
    value <- themode
    scalevalue <- therange[2] - therange[1]
    percentile.value <- therange[2]
    percentile.value <- percentile.value - 0.001
  }
  if (is.null(themean) && is.null(themode)) {
    stopifnot((themedian <= therange[2]) | (themedian >= therange[1]))
  } else if (is.null(themean) && is.null(themedian)) {
    stopifnot((themode <= therange[2]) | (themode >= therange[1]))
  }


  if (!is.null(themode) || !is.null(themedian)) {
    pr_n <- percentile
  }
  if (!is.null(themean)) {
    pr_n <- 0.999
  }




  if (is.null(themode) && is.null(themedian)) {
    to.minimize <- function(a) {
      abs(qbeta(pr_n, shape1 = a, shape2 = a * (1 - themean) / themean) - percentile.value)
    }
  } else if (is.null(themean) && is.null(themode)) {
    to.minimize <- function(a) {
      abs(qbeta(pr_n, shape1 = a, shape2 = (3 * a * (1 - themedian) + 2 * themedian - 1) / (3 * themedian)) - percentile.value)
    }
  } else if (is.null(themean) && is.null(themedian)) {
    to.minimize <- function(a) {
      abs(qbeta(pr_n, shape1 = a, shape2 = (a * (1 - themode) + 2 * themode - 1) / themode) - percentile.value)
    }
  }


  estimate <- optim(runif(1, 1, 10), to.minimize, lower = 0.1, upper = 10^4, method = "Brent")


  finalshape1 <- estimate$par

  if (is.null(themode) && is.null(themedian)) {
    finalshape2 <- finalshape1 * (1 - themean) / themean
  } else if (is.null(themean) && is.null(themode)) {
    finalshape2 <- (3 * finalshape1 * (1 - themedian) + 2 * themedian - 1) / (3 * themedian)
  } else if (is.null(themean) && is.null(themedian)) {
    finalshape2 <- (finalshape1 * (1 - themode) + 2 * themode - 1) / themode
  }

  set.seed(seed)
  sample_beta <- rbeta(nsims, finalshape1, finalshape2)
  param <- c(a = finalshape1, b = finalshape2)
  input <- c(tmetric = value, percentile = percentile, scalevalue=scalevalue, percentile.value = percentile.value)
  names(input)[1] <- name

  out <- list(parameters = param, summary = summary(sample_beta), input = input)
  print(out)
  class(out) <- "PriorGen"
  return(out)
}
