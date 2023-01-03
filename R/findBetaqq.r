#' The findbetaqq function
#'
#' A function to estimate the parameters alpha and beta of a Beta distribution based on the existing prior belief (data and/or expert opinion) about the values of two distinct percentiles.
#' 
#'
#' @usage findbeta(themean=NULL, themedian=NULL, themode=NULL, percentile=0.95,lower.v=F, percentile.value)
#'
#' @param percentile.value1: specify the value for the first percentile. It takes a value between 0 and 1.
#' @param percentile1: specify which is the percentile that corrersponds to percentile.value1. It takes a value between 0 and 1.
#' @param percentile.value2: specify the value for the second percentile. It takes a value between 0 and 1.
#' @param percentile2: specify which is the percentile that corrersponds to percentile.value2. It takes a value between 0 and 1.
#'
#' @examples
#' ## We believe that 20% of the units in an area/region have a prevalence of
#' ## disease/infection less than or equal to 0.30 while at the same time we are 90%
#' ## certain that the prevalence is less than 0.60
#'
#' findbetaqq(percentile.value1=0.30,percentile1=0.20,percentile.value2=0.60,percentile2=0.90)
#' 
#' @export 
#' @param parameters: The beta distribution parameters Beta(a,b)
#' @param summary: A basic summary of the elicited prior
#' @param input: The initial input value that produced the above prior.
#'
#' @references
#' Kostoulas, P., Nielsen, S. S., Branscum, A. J., Johnson, W. O., Dendukuri, N., Dhand, N. K., Toft, N., Gardner, I. A. (2017): Reporting guidelines for diagnostic accuracy studies that use Bayesian latent class models (STARD–BLCM). Statistics in medicine, 23, 3603–3604.


findbetaqq<-function(percentile.value1,percentile1,percentile.value2,percentile2, seed=280385){
  time=proc.time()
  
  findcentiles<-function(x)
  {
    c(
      F1=qbeta(percentile1, x[1], x[2]) - percentile.value1,
      F2=qbeta(percentile2, x[1], x[2]) - percentile.value2
    )
  }
  #library(rootSolve)
  
  k=1
  ss <- suppressWarnings(multiroot(f = findcentiles, start = c(2,2),maxiter = 1000))
       finalshape1=ss$root[1]
       finalshape2=ss$root[2]
       set.seed(seed)
       sample_beta=rbeta(10000,finalshape1,finalshape2)
  # while((ss$estim.precis=="NaN" | ss$estim.precis>1) | k>=100){
  #   if(ss$estim.precis!="NaN" | ss$estim.precis<=1){
  #     finalshape1=ss$root[1]
  #     finalshape2=ss$root[2]
  #     sample_beta=rbeta(10000,finalshape1,finalshape2)
  #   }else{
  #     ss <- multiroot(f = findcentiles, start = c(sample(1e+8,1),sample(1e+18,1)),maxiter = 1000)
  #     finalshape1=ss$root[1]
  #     finalshape2=ss$root[2]
  #     sample_beta=rbeta(10000,finalshape1,finalshape2)
  #   }
  #   k=k+1
  # }
       
  param=c(a=finalshape1,b=finalshape2)
  input=c(percentile.value1=percentile.value1,percentile1=percentile1,
          percentile.value2=percentile.value2,percentile2=percentile2)
  
  out<-list(parameters=param,summary=summary(sample_beta),input=input)
  class(out)<-"PriorGen"
  invisible(return(out))
}
