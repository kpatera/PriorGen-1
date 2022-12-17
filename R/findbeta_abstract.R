#' The findbeta (abstract) function
#'
#' A function to estimate the parameters alpha and beta of a Beta distribution based on the existing prior beliefs (data and/or expert opinion). 
#' General information is provided about the mean in terms of c("Very low","Low","Average","High","Very high"). The same holds for the variance parameter.
#'
#'
#' @usage function(themean.cat=c("Very low","Low","Average","High","Very high"),
#' thevariance.cat=c("Very low","Low","Average","High","Very high"), 
#' silent=TRUE, seed=280385, nsims=10000)
#' 
#' @param themean.cat: specify your prior belief about the mean. It takes a value among c("Very low","Low","Average","High","Very high").
#' @param thevariance.cat: specify your prior belief about the variance. It takes a value among c("Very low","Low","Average","High","Very high").
#' @param silent: If TRUE an extended output is printed. If FALSE and stored in an object the function runs silently.
#' @param seed: A fixed seed for replication purposes.
#' @param nsims: Number of simulations for the creation of various summary metrics of the elicited prior.
#'
#' @examples
#' ##Example 1
#' ##Based on the available literature the mean value for the sensitivity of a test
#' ##is expected to be generally low and its variance not that low but not that much neither.
#' 
#' findbeta_abstract(themean.cat="Low", thevariance.cat="Average")
#' 
#' 
#' @export 
#' @param parameters: The beta distribution parameters Beta(a,b)
#' @param summary: A basic summary of the elicited prior
#' @param input: The initial input value that produced the above prior.
#'
#' @references
#' Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2005): Estimation of diagnostic test sensitivity and specificity through Bayesian modeling. Preventive veterinary medicine, \bold{68}, 145--163.

findbeta_abstract<-function(themean.cat=c("Very low","Low","Average","High","Very high"),
                       thevariance.cat=c("Very low","Low","Average","High","Very high"),
                       silent=TRUE, seed=280385,nsims=10000){
  levels=c("Very low","Low","Average","High","Very high")
  themean=seq(0.1,0.9,length.out=5)[which(levels==themean.cat)]
  alpha=0.99995
  thevariance.optim=seq(0.001,0.5,by=0.001)
  thevar_max=thevariance.optim[length(which(themean+qnorm(alpha)*thevariance.optim<1))]
  thevariance=seq(0.001,thevar_max,length.out=5)[which(levels==thevariance.cat)]
  out1=NULL
  percentile=0.9999
  name<-"themean";value=themean;scalevalue=thevariance
  percentile.value=themean+qnorm(alpha)*thevariance

  
  
  a=runif(1,1,10)
  pr_n=0.999
  
  to.minimize<-function(a){
    abs(qbeta(pr_n, shape1=a, shape2=a * (1 - themean) / themean) - percentile.value)
    }

  estimate<-optim(runif(1,1,10),to.minimize,lower = 0.1, upper = 10^4, method = "Brent")
  
  
  finalshape1=estimate$par
  finalshape2 =finalshape1 * (1 - themean) / themean
  
  set.seed(seed)
  sample_beta=rbeta(nsims,finalshape1,finalshape2)
  param=c(a=finalshape1,b=finalshape2)
  input=c(tmetric=value,scalemetric=scalevalue,percentile.value=pr_n); names(input)[1]<-name
  names(input)[1]<-name
  
  if(silent==FALSE){
    cat(paste("The desired Beta distribution that satisfies the specified conditions is: Beta(", round(finalshape1,2), round(finalshape2,2),") \n",
              "Verification: The percentile value",round(qbeta(pr_n, finalshape1, finalshape2),2), "corresponds to the",pr_n,"th percentile \n"))
    cat(paste("Descriptive statistics for this distribution can be found below (or in the defined object):\n"))
    return(list(parameters=param,summary=summary(sample_beta),input=input))
  }
  invisible(return(list(parameters=param,summary=summary(sample_beta),input=input)))
  
}