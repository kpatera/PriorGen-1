#' The findbetamupsi (raw) function
#'
#' A function to estimate (a) the parameters of a Beta distribution for the expected mean of a proportion - usually the prevalence of disease/infection for the units in an area/region and (b) the parameters of a Gamma distribution expressing our prior belief about the variability of the prevalence estimates across the units of the area/region under consideration.
#'
#'
#' @usage findbetamupsi_raw(themean = 0.2, thevariance = 0.05, thepsi = 0.5,
#' seed = 280385, nsims = 10000, root.method = "multiroot")
#'
#' @param themean specify your prior belief about the mean. It takes a value between 0 and 1.
#' @param thevariance specify your prior belief about thevariance. If the selected variance is larger than possible, the variance will be adjusted downwards to create comply with the range of a probability.
#' @param thepsi specify your prior belief about the between-study variance. Large values of psi correspond to low heterogeneity, low values of psi correspond to large heterogeneity.
#' @param seed A fixed seed for replication purposes.
#' @param nsims Number of simulations for the creation of various summary metrics of the elicited prior.
#' @param root.method Choose between two alternatives to solve the two non-linear equations to identify the hyperparameters of psi. root.method="multiroot" involves the basic function of the rootSolve package, root.method="nleqslv" involves the base  functions of the nleqslv package.
#'
#' @examples
#' ## Example
#' ## The mean prevalence of a disease/infection for the units
#' ## within an area/region is thought to be 0.20 with a variance
#' ## equal to 0.03, with large heterogeneity i.e. psi equals to 0.15.
#'
#' findbetamupsi_raw(themean = 0.20, thevariance = 0.03, thepsi = 0.15)
#' findbetamupsi_raw(
#'   themean = 0.20, thevariance = 0.03, thepsi = 0.15,
#'   root.method = "nleqslv"
#' )
#'
#' @export
#' @return param_beta: The beta distribution parameters Beta(a,b)
#' @return param_gamma: The gamma distribution parameters gamma(a,b)
#' @return summary: A basic summary of the elicited prior
#' @return input: The initial input value that produced the above prior.
#' @return param_upper: simulated mu and psi of Beta(mu psi,psi(1-mu))
#'
#' @references
#' Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2005): Estimation of diagnostic test sensitivity and specificity through Bayesian modeling. Preventive veterinary medicine, \bold{68}, 145--163.

findbetamupsi_raw <- function(themean = 0.2, thevariance = 0.05, thepsi = 0.5,
                              seed = 280385, nsims = 10000, root.method = "multiroot") {
  # themean=0.2
  # thevariance = 0.05
  # psi.percentile=0.90
  # percentile.median=0.50
  # percentile95value=0.60
  # silent = TRUE
  alpha <- 0.9995
  pr_n <- 0.9999
  percentile.value <- themean + qnorm(alpha) * thevariance
  if ((themean + qnorm(alpha) * thevariance) > 1) {
    c("The set variance resulted in percentiles lower than 0 or greater than 1.
      The variance was adjusted.")
  }
  if ((themean + qnorm(alpha) * thevariance) > 1) {
    thevariance <- (1 - themean) / qnorm(alpha)
    percentile.value <- themean + qnorm(alpha) * thevariance
  }
  if ((themean - qnorm(alpha) * thevariance) < 0) {
    thevariance <- (themean) / (qnorm(alpha))
    percentile.value <- themean + qnorm(alpha) * thevariance
  }

  stopifnot(themean <= percentile.value)
  stopifnot(thepsi >= 0 && thepsi < 3)

  a <- runif(1, 1, 10)
  to.minimize <- function(a) {
    abs(qbeta(pr_n, shape1 = a, shape2 = a * (1 - themean) / themean) - percentile.value)
  }

  estimate <- optim(runif(1, 1, 10), to.minimize, lower = 0.1, upper = 10^4, method = "Brent")
  alpha_mu <- estimate$par
  beta_mu <- alpha_mu * (1 - themean) / themean
  mu <- alpha_mu / (alpha_mu + beta_mu)
  # gu <- alpha_mu + beta_mu

  alpha_t <- mu * thepsi
  beta_t <- thepsi * (1 - mu)

  model <- function(x) {
    c(F1 = qgamma(0.5, shape = x[1], scale = 1 / x[2]) -
      alpha_t, F2 = qgamma(0.05, shape = x[1], scale = 1 / x[2]) -
      beta_t)
  }

  if (root.method == "multiroot") {
    ss2 <- multiroot(f = model, start = c(5, 2), positive = T)
    a <- rgamma(nsims, ss2$root[1], ss2$root[2])
  } else if (root.method == "nleqslv") {
    ss2 <- nleqslv(c(5, 2), model, control = list(btol = .01))
    a <- rgamma(nsims, ss2$x[1], ss2$x[2])
  }

  b <- rbeta(nsims, alpha_mu, beta_mu)
  param <- c(a = alpha_mu, b = beta_mu)
  sample_beta <- rbeta(nsims, a * b, a * (1 - b))

  input <- c(
    themean = themean, thevariances = thevariance, percentile = pr_n,
    percentile.value = percentile.value, thepsi = thepsi
  )

  out <- list(param_beta = param, param_gamma = ss2$root, summary = summary(sample_beta), input = input, param_upper = list(at = a * b, bt = a * (1 - b)))
  class(out) <- "PriorGen"
  invisible(return(out))
}
