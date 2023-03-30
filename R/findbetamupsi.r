#' The findbetamupsi function
#'
#' A function to estimate (a) the parameters of a Beta distribution for the expected mean of a proportion - usually the prevalence of disease/infection for the units in an area/region and (b) the parameters of a Gamma distribution expressing our prior belief about the variability of the prevalence estimates across the units of the area/region under consideration.
#'
#'
#' @usage findbetamupsi(themean, percentile=0.95, lower.v=T, percentile.value,  psi.percentile=0.90, percentile.median, percentile95value)
#'
#' @param themean: specify your prior belief about the mean. It takes a value between 0 and 1.
#' @param percentile: specify the level of confidence that the true value of the mean is greater or lower than the percentile.value. It takes a value between 0 and 1 and the default is 0.95.
#' @param lower.v: logical, if TRUE the specified percentile.value is the upper limit for the mean at the specified confidence level (percentile). If FALSE the specified percentile.value is the lower limit for the mean at the specified confidence level (percentile).The default is TRUE.
#' @param percentile.value: specify the upper or lower limit for the mean at the specified level of confidence (percentile). It takes a value between 0 and 1.
#' @param psi.percentile: specify the level of confidence that a certain fraction of the units under study has a prevalence less than the percentile.median. It takes a value between 0 and 1 and the default is 0.90.
#' @param percentile.median: specify the median value that corresponds to the defined psi.percentile. It takes a value between 0 and 1 and has to be higher than both themean and the percentile.
#' @param percentile95value: specify the value that the percentile.median does not exceed with 95% confidence. It takes a value between 0 and 1 and has to be higher than the percentile.median.
#' @param seed: A fixed seed for replication purposes.
#' @param nsims: Number of simulations for the creation of various summary metrics of the elicited prior.
#' @param root.method: Choose between two alternatives to solve the two non-linear equations to identify the hyperparameters of psi. root.method="multiroot" involves the basic function of the rootSolve package, root.method="nleqslv" involves the base functions of the nleqslv package.
#'
#' @examples
#' ## Example
#' ## The mean prevalence of a disease/infection for the units within an area/region
#' ## is thought to be 0.20 and we are 99% confident that it is not more than 0.40.
#' ## Within this area, we are also confident that 90% of all units have a prevalence
#' ## less or equal to 0.50 and we are 95% certain that it does not exceed 0.60
#'
#' res1 <- findbetamupsi(themean = 0.20, percentile = 0.99, lower.v = TRUE, percentile.value = 0.30, psi.percentile = 0.90, percentile.median = 0.50, percentile95value = 0.60)
#' res2 <- findbetamupsi(themean = 0.20, percentile = 0.99, lower.v = TRUE, percentile.value = 0.30, psi.percentile = 0.90, percentile.median = 0.50, percentile95value = 0.60, root.method = "nleqslv")
#' res1
#' res2
#'
#' @export
#' @param param_beta: The beta distribution parameters Beta(a,b)
#' @param param_gamma: The gamma distribution parameters gamma(a,b)
#' @param summary: A basic summary of the elicited prior
#' @param input: The initial input value that produced the above prior.
#' @param param_upper: simulated mu and psi of Beta(mu psi,psi(1-mu))
#'
#' @references
#' Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2005): Estimation of diagnostic test sensitivity and specificity through Bayesian modeling. Preventive veterinary medicine, \bold{68}, 145--163.

findbetamupsi <- function(themean, percentile = 0.95, lower.v = T, percentile.value,
                          psi.percentile = 0.90, percentile.median, percentile95value,
                          seed = 280385, nsims = 10000, root.method = "multiroot") {
  if (lower.v == T) {
    pr_n <- percentile
  } else {
    pr_n <- 1 - percentile
  }
  stopifnot((lower.v == T && themean <= percentile.value) |
    (lower.v == F && themean >= percentile.value))
  stopifnot(percentile.median > themean && percentile95value >
    percentile.median)
  stopifnot((lower.v == T && percentile.value < percentile.median) |
    (lower.v == F && percentile.value != percentile.median))
  stopifnot(psi.percentile > 0.5 && psi.percentile < 1)
  a <- runif(1, 1, 10)
  to.minimize <- function(a) {
    abs(qbeta(pr_n, shape1 = a, shape2 = a * (1 - themean) / themean) -
      percentile.value)
  }
  estimate <- optim(runif(1, 1, 10), to.minimize,
    lower = 0.1,
    upper = 10^4, method = "Brent"
  )
  alpha_mu <- estimate$par
  beta_mu <- alpha_mu * (1 - themean) / themean
  mu <- alpha_mu / (alpha_mu + beta_mu)
  gu <- alpha_mu + beta_mu
  f <- function(ga) {
    abs(qbeta(psi.percentile, mu * ga, ga *
      (1 - mu)) - percentile.median)
  }
  r <- optim(1, f, lower = 0, upper = 10^4, method = "Brent")
  alpha_t <- r$par[1]
  f <- function(ga) {
    abs(qbeta(psi.percentile, mu * ga, ga *
      (1 - mu)) - percentile95value)
  }
  r <- optim(1, f, lower = 0, upper = 10^4, method = "Brent")
  beta_t <- r$par[1]
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
    themean = themean, percentile = percentile,
    percentile.value = percentile.value, psi.percentile = psi.percentile,
    percentile.median = percentile.median, percentile95value = percentile95value
  )

  out <- list(param_beta = param, param_gamma = ss2$root, summary = summary(sample_beta), input = input, param_upper = list(at = a * b, bt = a * (1 - b)))
  class(out) <- "PriorGen"
  invisible(return(out))
}
