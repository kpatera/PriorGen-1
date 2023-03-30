##' #' The findbeta lines function
#'
#' A function that plots any object on top of another plot of the class findbeta.
#'
#' @usage lines(findbeta.object, ...)
#' @rdname lines.PriorGen
#' @param findbeta.object: An object of type findbeta produces of one of the other PriorGen functions.
#'
#' @examples
#' ##Example 1
#' ##Based on the available literature the mean value for the sensitivity of a test
#' ##is expected to be generally low and its variance not that low but not that much neither.
#' 
#' res_abs_1<-findbeta_abstract(themean.cat="Low", thevariance.cat="Average")
#' plot(res_abs_1,main="Plot of the findbeta_abstract function",lwd=3,ylim=c(0,7))
#' res_abs_2<-findbeta_abstract(themean.cat="High", thevariance.cat="Average")
#' lines(res_abs_2,lwd=3)
#' 
#' ## Example 2 
#' ## Hierarchical prior
#' 
#' res_mult_1=findbetamupsi(themean=0.10, percentile=0.79, lower.v=TRUE, percentile.value=0.26, psi.percentile=0.95, percentile.median=0.28, percentile95value=0.3)
#' plot(res_mult_1,main="Plot of the findbetamupsi function",lwd=3,ylim=c(0,7))
#' res_mult_2=findbetamupsi(themean=0.12, percentile=0.79, lower.v=TRUE, percentile.value=0.26, psi.percentile=0.95, percentile.median=0.28, percentile95value=0.3)
#' lines(res_mult_2,lwd=3)
#' @export 


lines.PriorGen<-function(findbeta.object,...){
  if(class(findbeta.object)!="PriorGen")
    stop("Provide an object of class PriorGen")
  
  if(length(findbeta.object)==3){
    a_plot=findbeta.object$parameters[1]
    b_plot=findbeta.object$parameters[2]
    x=seq(0,1,0.001)
    y=dbeta(seq(0,1,0.001),a_plot,b_plot)
    lines(x,y,...)
  }else{
    a_plot=findbeta.object$parameters$at
    b_plot=findbeta.object$parameters$bt
    x=seq(0,1,0.001)
    y=dbeta(seq(0,1,0.001),a_plot,b_plot)
    lines(x,y,...)
  }
}