#' The findbeta (raw) function
#'
#' A function to estimate the parameters alpha and beta (a,b) of a Beta distribution based on the existing prior beliefs (data and/or expert opinion). 
#' Information should be provided on the raw values of the mean (or the median or the mode) and a corresponding scale metric, either the variance or the range of the parameter.
#'
#' @usage findbeta_raw(themean=NULL,themedian=NULL,themode=NULL,
#'  thevariance=NULL, therange=c(0,1), seed=280385, nsims=10000)
#'
#' @param themean specify your prior belief about the mean. It takes a value between 0 and 1. Not to be specified if a value has been given for the median or the mode.
#' @param themedian specify your prior belief about the median. It takes a value between 0 and 1. Not to be specified if a value has been given for the mean or the mode.
#' @param themode specify your prior belief about the mode. It takes a value between 0 and 1. Not to be specified if a value has been given for the mean or the median.
#' @param thevariance specify your prior belief about the variance. If the selected variance is larger than possible, the variance will be adjusted downwards to create comply with the range of a probability.
#' @param therange specify your prior belief about the range. It should be a two number vector that c(ul,ll), where ul>0, ll<1 and ul<ll. This scale metric applies for themode and themedian options.
#' @param seed A fixed seed for replication purposes.
#' @param nsims Number of simulations for the creation of various summary metrics of the elicited prior.
#'
#' @examples
#' ## Example 1
#' ## Based on the available literature the mean value for the sensitivity of a test
#' ## is expected to be 0.90 and its variance equals to 0.1.
#'
#' findbeta_raw(themean = 0.90, thevariance = 0.1)
#'
#' ## Example 2
#' ## Based on the available literature the median value for the specificity of a
#' ## test is expected to be 0.99 and its range between 0.1 and 1.
#'
#' findbeta_raw(themedian = 0.70, therange = c(0.1, 1))
#'
#'
#' # Mode
#'
#' findbeta_raw(themode = 0.70, therange = c(0.1, 1))
#'
#' @export
#' @return parameters The beta distribution parameters Beta(a,b)
#' @return summary A basic summary of the elicited prior
#' @return input The initial input value that produced the above prior.
#'
#' @references
#' Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2005): Estimation of diagnostic test sensitivity and specificity through Bayesian modeling. Preventive veterinary medicine, \bold{68}, 145--163.


findbeta_raw <- function(themean = NULL, themedian = NULL, themode = NULL,
                         thevariance = NULL, therange = c(0, 1), seed = 280385,
                         nsims = 10000) {
  # Cannot handle Beta(a,b) with a,b<1
  
  if((is.null(themean) + is.null(themedian) + is.null(themode)) == 2) 
    stop("Error: Input at least one of the following: 1. themean, 2. themedian, 3. themode.")
  
  if(((thevariance <= 0.5) + (therange[2] <= 1) + (therange[1] >= 0)) == 3) 
    stop("Error: Input at least one of the following: 1. thevariance (for themean), 2. therange (for themedian/themode)")

  alpha <- 0.9999995
  pr_n <- 0.999
  out <- NULL

  if (!is.null(themean)) {
    if ((themean + qnorm(alpha) * thevariance) > 1) {
      c("The set variance resulted in percentiles lower than 0 or greater than 1. The variance was adjusted.")
    }
    if ((themean + qnorm(alpha) * thevariance) > 1) {
      thevariance <- (1 - themean) / qnorm(alpha)
      percentile.value <- themean + qnorm(alpha) * thevariance
    }
    if ((themean - qnorm(alpha) * thevariance) < 0) {
      thevariance <- (themean) / (qnorm(alpha))
      percentile.value <- themean + qnorm(alpha) * thevariance
    }
  }

  if (!is.null(themean)) {
    name <- "themean"
    value <- themean
    scalevalue <- thevariance
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


  estimate <- optim(runif(1, 1, 10), to.minimize, lower = 0.001, upper = 10^4, method = "Brent")


  a_p <- estimate$par

  if (is.null(themode) && is.null(themedian)) {
    b_p <- a_p * (1 - themean) / themean
  } else if (is.null(themean) && is.null(themode)) {
    b_p <- (3 * a_p * (1 - themedian) + 2 * themedian - 1) / (3 * themedian)
  } else {
    b_p <- (a_p * (1 - themode) + 2 * themode - 1) / themode
  }

  set.seed(seed)
  sample_beta <- rbeta(nsims, a_p, b_p)
  param <- c(a = a_p, b = b_p)
  input <- c(tmetric = value, scalemetric_var_or_range = scalevalue)
  names(input)[1] <- name
  names(input)[1] <- name

  out <- list(parameters = param, summary = summary(sample_beta), input = input)
  class(out) <- "PriorGen"
  return(out)
}
