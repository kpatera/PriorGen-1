# The findbeta plot function
#
# A function that plots any object of the class findbeta.
#
# @usage function(findbeta.object, lines=FALSE, ...)
# 
# @param findbeta.object: An object of type findbeta produces of one of the other PriorGen functions.
# @param lines: If lines=TRUE then it plots over the last plotted plot.
# @param ...: Additional base plot parameters.
#
# @examples
# ##Example 1
# ##Based on the available literature the mean value for the sensitivity of a test
# ##is expected to be generally low and its variance not that low but not that much neither.
# 
# res_abs_1<-findbeta_abstract(themean.cat="Low", thevariance.cat="Average")
# findbeta_plot(res_abs_1,main="Plot of the findbeta_abstract function",lwd=3)
# 
# @export 
# @param parameters: The beta distribution parameters Beta(a,b)
# @param summary: A basic summary of the elicited prior
# @param input: The initial input value that produced the above prior.
#
# @references
# Branscum, A. J., Gardner, I. A., & Johnson, W. O. (2005): Estimation of diagnostic test sensitivity and specificity through Bayesian modeling. Preventive veterinary medicine, \bold{68}, 145--163.
plot.PriorGen<-function(findbeta.object,...){
  a_plot=findbeta.object$parameters[1]
  b_plot=findbeta.object$parameters[2]
  x=seq(0,1,0.001)
  y=dbeta(seq(0,1,0.001),a_plot,b_plot)
  plot(x,y,...)
#  if(lines==F){plot(x,x_beta,...)
#  }else if (lines==T){lines(x,x_beta,...)
 # }
}

lines.PriorGen<-function(findbeta.object,...){
  a_plot=findbeta.object$parameters[1]
  b_plot=findbeta.object$parameters[2]
  x=seq(0,1,0.001)
  y=dbeta(seq(0,1,0.001),a_plot,b_plot)
  lines(x,y,...)
  #  if(lines==F){plot(x,x_beta,...)
  #  }else if (lines==T){lines(x,x_beta,...)
  # }
}



plot.PriorGen2<-function(findbeta.object,...){
  a_plot=mean(findbeta.object$param_upper[1])
  b_plot=mean(findbeta.object$param_upper[2])
  x=seq(0,1,0.001)
  y=dbeta(seq(0,1,0.001),a_plot,b_plot)
  plot(x,y,...)
  #  if(lines==F){plot(x,x_beta,...)
  #  }else if (lines==T){lines(x,x_beta,...)
  # }
}

lines.PriorGen2<-function(findbeta.object,...){
  a_plot=findbeta.object$parameters[1]
  b_plot=findbeta.object$parameters[2]
  x=seq(0,1,0.001)
  y=dbeta(seq(0,1,0.001),a_plot,b_plot)
  lines(x,y,...)
  #  if(lines==F){plot(x,x_beta,...)
  #  }else if (lines==T){lines(x,x_beta,...)
  # }
}

print.PriorGen<-function(findbeta.object,...){
  cat(paste("The desired Beta distribution that satisfies the specified conditions is: Beta(", round(findbeta.object$parameters[1],2),",", round(findbeta.object$parameters[2],2),") \n",
            "Verification: The percentile value",round(qbeta(findbeta.object$input[2], findbeta.object$parameters[1], findbeta.object$parameters[2]),2), "corresponds to the",findbeta.object$input[2],"th percentile \n"))
  cat(paste("Descriptive statistics for this distribution can be found below:\n"))
  print(findbeta.object)
}

print.PriorGen2<-function(findbeta.object,...){
  print (paste("The desired Beta distribution that satisfies the specified conditions about the mean of the prevalence 'mu' is: Beta(", round(findbeta.object$param_beta[1],2), round(findbeta.object$param_beta[2],2),")"))
  print (paste("The desired Gamma distribution that satisfies the specified conditions about the variability 'psi' of the prevalence is: Gamma(", round(findbeta.object$param_gamma$root[1],2), round(findbeta.object$param_gamma$root[2],2),")"))
  print("Descriptive statistics for this distrubiton are:")
  print(findbeta.object)  
}