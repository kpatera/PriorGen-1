##' #' The findbeta plot function
#'
#' A function that prints a summary any object of the class findbeta.
#'
#' @usage print(findbeta.object)
#' @rdname print.PriorGen
#' @param findbeta.object: An object of type findbeta produced of one of the main PriorGen functions.
#'
#' @examples
#' ##Example 1
#' ##Based on the available literature the mean value for the sensitivity of a test
#' ##is expected to be generally low and its variance not that low but not that much neither.
#' res_abs_1<-findbeta_abstract(themean.cat="Low", thevariance.cat="Average")
#' print(res_abs_1)
#' 
#' Example 2 
#' Hierarchical prior
#' res_mult_1=findbetamupsi(themean=0.10, percentile=0.79, lower.v=TRUE, percentile.value=0.26, psi.percentile=0.95, percentile.median=0.28, percentile95value=0.3)
#' print(res_mult_1)
#' 
#' @export
print.PriorGen<-function(findbeta.object,...){
  if(class(findbeta.object)!="PriorGen")
    stop("Provide an object of class PriorGen")
  if(length(findbeta.object)==3){
    cat(paste("The desired Beta distribution that satisfies the specified conditions is: Beta(", round(findbeta.object$parameters[1],2),",", round(findbeta.object$parameters[2],2),") \n",
              "Verification: The percentile value",round(qbeta(findbeta.object$input[2], findbeta.object$parameters[1], findbeta.object$parameters[2]),2), "corresponds to the",findbeta.object$input[2],"th percentile \n"))
    cat(paste("Descriptive statistics for this distribution can be found below:\n"))
  }else{
    print (paste("The desired Beta distribution that satisfies the specified conditions about the mean of the prevalence 'mu' is: Beta(", round(findbeta.object$param_beta[1],2), round(findbeta.object$param_beta[2],2),")"))
    print (paste("The desired Gamma distribution that satisfies the specified conditions about the variability 'psi' of the prevalence is: Gamma(", round(findbeta.object$param_gamma$root[1],2), round(findbeta.object$param_gamma$root[2],2),")"))
    print("Descriptive statistics for this distrubiton are:")
    
  }
  print(findbeta.object)
}