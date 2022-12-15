findbeta_abstract<-function(themean.cat=c("Very low","Low","Average","High","Very high"),
                       thevariance.cat=c("Very low","Low","Average","High","Very high"),
                       dens=FALSE){
  levels=c("Very low","Low","Average","High","Very high")
  themean=seq(0,1,length.out=5)[which(levels==themean.cat)]
  alpha=0.99995
  thevariance.optim=seq(0,0.5,by=0.001)
  thevar_max=thevariance.optim[length(which(themean+qnorm(alpha)*thevariance.optim<1))]
  thevariance=seq(0,thevar_max,length.out=5)[which(levels==thevariance.cat)]
  
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
  
  set.seed(280385)
  sample_beta=rbeta(10000,finalshape1,finalshape2)
  
  cat(paste("The desired Beta distribution that satisfies the specified conditions is: Beta(", round(finalshape1,2), round(finalshape2,2),") \n",
            "Verification: The percentile value",round(qbeta(pr_n, finalshape1, finalshape2),2), "corresponds to the",pr_n,"th percentile \n"))
  cat(paste("Descriptive statistics for this distribution can be found below (or in the defined object):\n"))
  
  param=c(a=finalshape1,b=finalshape2)
  input=c(tmetric=value,scalemetric=scalevalue,percentile.value=pr_n); names(input)[1]<-name
  x_beta=dbeta(seq(0,1,0.001),finalshape1,finalshape2)
  if(dens==FALSE){density_plot=NA}else if(dens==TRUE){density_plot=plot(x,x_beta,type = "l",lwd=3, 
                                                                        xlab="Elicited beta prior",
                                                                        ylab = "Density")
  }
  return(list(parameters=param,summary=summary(sample_beta),input=input,figure=density_plot))
  
}


findbeta_abstract(themean = 0.5,thevariance = 0.2,dens = T) # Check algorithm for limit of variance!!
temp=findbeta_abstract(themean = 0.8,thevariance = 0.2,dens = F)
temp$parameters
temp$summary
temp$input

