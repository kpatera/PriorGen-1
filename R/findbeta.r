#' The findbeta function
#'
#' A function to estimate the parameters alpha and beta of a Beta distribution based on the existing prior beliefs (data and/or expert opinion). 
#' Information should be provided about the mean (or the median or the mode) and whether it is lower or greater that a certain value with a pre-specified certainty (usually set at 95\%)
#'
#' @usage findbeta(themean=NULL, themedian=NULL, themode=NULL,
#'  percentile=0.95,lower.v=F, percentile.value,
#'  seed=280385, nsims=10000)
#'
#' @param themean specify your prior belief about the mean. It takes a value between 0 and 1. Not to be specified if a value has been given for the median or the mode.
#' @param themedian specify your prior belief about the median. It takes a value between 0 and 1. Not to be specified if a value has been given for the mean or the mode.
#' @param themode specify your prior belief about the mode. It takes a value between 0 and 1. Not to be specified if a value has been given for the mean or the median.
#' @param percentile specify the level of confidence that the true value of the mean (or the median or the mode) is greater or lower than the percentile.value. It takes a value between 0 and 1 and the default =0.95.
#' @param lower.v logical, if TRUE the specified percentile.value is the upper limit for the mean (or the median or the mode) at the specified confidence level (percentile). If {FALSE} the specified percentile.value is the lower limit for the mean (or the median or the mode) at the specified confidence level (percentile). The default is FALSE.
#' @param percentile.value specify the upper or lower limit for the mean (or the median or the mode) at the specified level of confidence (percentile). It takes a value between 0 and 1.
#' @param seed A fixed seed for replication purposes.
#' @param nsims Number of simulations for the creation of various summary metrics of the elicited prior.
#'
#' @examples
#' ## Example 1
#' ## Based on the available literature the mean value for the sensitivity of a test
#' ## is expected to be 0.90 and we can be 95\% sure that it is higher than 0.80.
#'
#' findbeta(
#'   themean = 0.90, percentile = 0.95, lower.v = FALSE,
#'   percentile.value = 0.80, seed = 280385, nsims = 10000
#' )
#'
#' ## Example 2
#' ## Based on the available literature the median value for the specificity of a
#' ## test is expected to be 0.99 and we can be 95\% sure that it is higher than 0.90.
#'
#'findbeta(
#'  themedian = 0.99, percentile = 0.95, lower.v = FALSE,
#'  percentile.value = 0.90
#')
#'
#' ##Example 3
#' ##The most probable value (mode) for the prevalence of a disease/infection in a
#' ##population is expected to be 0.15 and we are 90\% sure that it is less than 0.40.
#'
#' findbeta(themode=0.15, percentile=0.90,lower.v=TRUE,
#' percentile.value=0.40)
#'
#' @return parameters: The beta distribution parameters Beta(a,b)
#' @return summary: A basic summary of the elicited prior
#' @return input: The initial input value that produced the above prior.
#'
#' @export
#' @references
#' Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2005): Estimation
#' of diagnostic test sensitivity and specificity through Bayesian modeling.
#' Preventive veterinary medicine, \bold{68}, 145--163.


findbeta <- function(themean = NULL, themedian = NULL, themode = NULL,
                     percentile = 0.95, lower.v = F, percentile.value,
                     seed = 280385, nsims = 10000) {
  stopifnot((is.null(themean) + is.null(themedian) + is.null(themode)) == 2)
  out <- NULL
  if (!is.null(themean)) {
    name <- "themean"
    value <- themean
  }
  if (!is.null(themedian)) {
    name <- "themedian"
    value <- themedian
  }
  if (!is.null(themode)) {
    name <- "themode"
    value <- themode
  }
  if (is.null(themode) && is.null(themedian)) {
    stopifnot((lower.v == T && themean <= percentile.value) | (lower.v == F && themean >= percentile.value))
  } else if (is.null(themean) && is.null(themode)) {
    stopifnot((lower.v == T && themedian <= percentile.value) | (lower.v == F && themedian >= percentile.value))
  } else {
    stopifnot((lower.v == T && themode <= percentile.value) | (lower.v == F && themode >= percentile.value))
  }

  if (lower.v == T) {
    pr_n <- percentile
  } else {
    pr_n <- 1 - percentile
  }

  if (is.null(themode) && is.null(themedian)) {
    to.minimize <- function(a) {
      abs(qbeta(pr_n, shape1 = a, shape2 = a * (1 - themean) / themean) - percentile.value)
    }
  } else if (is.null(themean) && is.null(themode)) {
    to.minimize <- function(a) {
      abs(qbeta(pr_n, shape1 = a, shape2 = (3 * a * (1 - themedian) + 2 * themedian - 1) / (3 * themedian)) - percentile.value)
    }
  } else {
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
  } else {
    finalshape2 <- (finalshape1 * (1 - themode) + 2 * themode - 1) / themode
  }

  set.seed(seed)
  sample_beta <- rbeta(nsims, finalshape1, finalshape2)
  param <- c(a = finalshape1, b = finalshape2)
  input <- c(tmetric = value, percentile = percentile, percentile.value = percentile.value)
  names(input)[1] <- name

  out <- list(parameters = param, summary = summary(sample_beta), input = input)
  class(out) <- c("PriorGen")
  return(out)
}

