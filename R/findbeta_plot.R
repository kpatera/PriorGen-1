#' The findbeta plot function
#'
#' A function that plots any object of the class findbeta.
#'
#' @usage function(findbeta.object, lines=FALSE, ...)
#' 
#' @param findbeta.object: An object of type findbeta produces of one of the other PriorGen functions.
#' @param lines: If lines=TRUE then it plots over the last plotted plot.
#' @param ...: Additional base plot parameters.
#'
#' @examples
#' ##Example 1
#' ##Based on the available literature the mean value for the sensitivity of a test
#' ##is expected to be generally low and its variance not that low but not that much neither.
#' 
#' res_abs_1<-findbeta_abstract(themean.cat="Low", thevariance.cat="Average")
#' findbeta_plot(res_abs_1,main="Plot of the findbeta_abstract function",lwd=3)
#' 
#' @export 
#' @param parameters: The beta distribution parameters Beta(a,b)
#' @param summary: A basic summary of the elicited prior
#' @param input: The initial input value that produced the above prior.
#' @import 
#'
#' @references
#' Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2005): Estimation of diagnostic test sensitivity and specificity through Bayesian modeling. Preventive veterinary medicine, \bold{68}, 145--163.
findbeta_plot<-function(findbeta.object,lines=FALSE,...){
  a_plot=findbeta.object$parameters[1]
  b_plot=findbeta.object$parameters[2]
  x=seq(0,1,0.001)
  x_beta=dbeta(seq(0,1,0.001),a_plot,b_plot)
  if(lines==F){plot(x,x_beta,...)
  }else if (lines==T){lines(x,x_beta,...)
  }
}
