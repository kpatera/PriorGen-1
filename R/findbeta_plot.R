findbeta_plot<-function(findbeta.object,lines=FALSE,...){
  finalshape1=findbeta.object$parameters[1]
  finalshape2=findbeta.object$parameters[2]
  x=seq(0,1,0.001)
  x_beta=dbeta(seq(0,1,0.001),finalshape1,finalshape2)
  if(lines==F){plot(x,x_beta,...)
  }else if (lines==T){lines(x,x_beta,...)
  }
}

# source("~/GitHub/PriorGen-1/R/findbeta_plot.r")
# fb1=findbeta(themean = 0.8,lower.v = T,percentile = 0.90,percentile.value = 0.95,silent = F)
# findbeta_plot(fb1,xlab="Elicited beta prior",ylab = "Density",lwd=3,type="l")
