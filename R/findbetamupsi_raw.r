#' The findbetamupsi (raw) function
#'
#' A function to estimate (a) the parameters of a Beta distribution for the expected mean of a proportion - usually the prevalence of disease/infection for the units in an area/region and (b) the parameters of a Gamma distribution expressing our prior beleif about the variability of the prevalence estimates across the units of the area/region under consideration.
#'
#'
#' @usage findbetamupsi_raw(themean, thevariance,  psi.percentile=0.90, percentile.median, percentile95value)
#' 
#' @param themean: specify your prior belief about the mean. It takes a value between 0 and 1.
#' @param thevariance: specify your prior belief about the variance. If the selected variance is larger than possible, the variance will be adjusted downwards to create comply with the range of a probability.
#' @param psi.percentile: specify the level of confidence that a certain fraction of the units under study has a prevalence less than the percentile.median. It takes a value between 0 and 1 and the default is 0.90.
#' @param percentile.median: specify the median value that corresponds to the defined psi.percentile. It takes a value between 0 and 1 and has to be higher than both themean and the percentile.
#' @param percentile95value: specify the value that the percentile.median does not exceed with 95% confidence. It takes a value between 0 and 1 and has to be higher than the percentile.median.
#'
#' @examples
#' ## Example
#' ## The mean prevalence of a disease/infection for the units within an area/region
#' ## is thought to be 0.20 with a variance equal to 0.03, 
#' ## we are also confident that 90% of all units have a prevalence
#' ## less or equal to 0.50 and we are 95% certain that it does not exceed 0.60
#'
#' findbetamupsi_raw(themean=0.20, thevariance = 0.03, psi.percentile=0.90, percentile.median=0.50, percentile95value=0.60)
#' 
#' @export 
#' @param parameters: The beta distribution parameters Beta(a,b)
#' @param bot_param: simulated mu and psi of Beta(mu psi,psi(1-mu))
#' @param summary: A basic summary of the elicited prior
#' @param input: The initial input value that produced the above prior.
#' 
#' @import 
#'
#' @references
#' Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2005): Estimation of diagnostic test sensitivity and specificity through Bayesian modeling. Preventive veterinary medicine, \bold{68}, 145--163.

findbetamupsi_raw<-function(themean=0.2, thevariance=0.05, 
                            psi.percentile=0.90, percentile.median, percentile95value,
                            seed=280385,silent=TRUE){
  # themean=0.2
  # thevariance = 0.05
  # psi.percentile=0.90
  # percentile.median=0.50
  # percentile95value=0.60
  # silent = TRUE
  alpha=0.9995
  pr_n = 0.9999
  percentile.value=themean+qnorm(alpha)*thevariance
  if((themean+qnorm(alpha)*thevariance)>1){out1=c("The set variance resulted in percentiles lower than 0 or greater than 1. The variance was adjusted.")}
  if((themean+qnorm(alpha)*thevariance)>1){thevariance=(1-themean)/qnorm(alpha);percentile.value=themean+qnorm(alpha)*thevariance}
  if((themean-qnorm(alpha)*thevariance)<0){thevariance=(themean)/(qnorm(alpha));percentile.value=themean+qnorm(alpha)*thevariance}
  
  stopifnot(themean <= percentile.value)
  stopifnot(psi.percentile > 0.5 && psi.percentile < 1)
  stopifnot(percentile.median > themean && percentile95value > percentile.median)
  stopifnot(percentile.value < percentile.median)
  
  a = runif(1, 1, 10)
  to.minimize <- function(a) {
    abs(qbeta(pr_n, shape1 = a, shape2 = a * (1 - themean)/themean) - percentile.value)
  }
  estimate <- optim(runif(1, 1, 10), to.minimize, lower = 0.1, upper = 10^4, method = "Brent")
  alpha_mu = estimate$par
  beta_mu = alpha_mu * (1 - themean)/themean
  mu = alpha_mu/(alpha_mu + beta_mu)
  gu = alpha_mu + beta_mu
  f <- function(ga) abs(qbeta(psi.percentile, mu * ga, ga * (1 - mu)) - percentile.median)
  r <- optim(1, f, lower = 0, upper = 10^4, method = "Brent")
  alpha_t = r$par[1]
  f <- function(ga) abs(qbeta(psi.percentile, mu * ga, ga * (1 - mu)) - percentile95value)
  r <- optim(1, f, lower = 0, upper = 10^4, method = "Brent")
  beta_t = r$par[1]
  model <- function(x) c(F1 = qgamma(0.5, shape = x[1], scale = 1/x[2]) - 
                           alpha_t, F2 = qgamma(0.05, shape = x[1], scale = 1/x[2]) - 
                           beta_t)
  ss2 <- multiroot(f = model, start = c(5, 2), positive = T)
  qgamma(0.5, ss2$root[1], ss2$root[2])
  qgamma(0.95, ss2$root[1], ss2$root[2])
  
  
  a=rgamma(10000,ss2$root[1],ss2$root[2])
  b=rbeta(10000,alpha_mu,beta_mu)
  param=c(a=alpha_mu,b=beta_mu)
  sample_beta=rbeta(10000,a*b,a*(1-b))
  
  input=c(themean=themean, thevariances=thevariance,percentile=pr_n,
          percentile.value=percentile.value, psi.percentile=psi.percentile, 
          percentile.median=percentile.median, percentile95value=percentile95value)
  
  if(silent==FALSE){
    print (paste("The desired Beta distribution that satisfies the specified conditions about the mean of the prevalence 'mu' is: Beta(", round(finalshape1,2), round(finalshape2,2),")"))
    print (paste("The desired Gamma distribution that satisfies the specified conditions about the variability 'psi' of the prevalence is: Gamma(", round(ss2$root[1],2), round(ss2$root[2],2),")"))
    #print ("The plot gives the specified prior beleif on the prevalence distribution.")
    #plot(density(rbeta(10000,a*b,a*(1-b))))
    print("Descriptive statistics for this distrubiton are:")
    return(list(parameters=param,bot_param=list(at=a*b,bt=a*(1-b)),summary=summary(sample_beta),input=input))
  }
  invisible(return(list(parameters=param,bot_param=list(at=a*b,bt=a*(1-b)),summary=summary(sample_beta),input=input)))
  
}



# temp<-findbetamupsi_raw(themean=0.30, thevariance=0.05, psi.percentile=0.90,
#                         percentile.median=0.6, percentile95value=0.70)
# temp$parameters
