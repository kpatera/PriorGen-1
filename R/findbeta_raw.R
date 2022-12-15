findbeta_raw<-function(themean=NULL,themedian=NULL,themode=NULL,
                   thevariance=NULL, therange=c(0,1), silent=TRUE, seed=280385,
                   nsims=10000){
  # Cannot handle Beta(a,b) with a,b<1
  stopifnot ( (is.null(themean) + is.null(themedian) + is.null(themode)) == 2)
  stopifnot ( ((thevariance<=0.5) + (therange[2]<=1) + (therange[1]>=0)) == 3)
  
  alpha=0.99995
  pr_n=0.9999
  out1=NULL
  
  if(!is.null(themean)){
    if((themean+qnorm(alpha)*thevariance)>1){out1=c("The set variance resulted in percentiles lower than 0 or greater than 1. The variance was adjusted.")}
    if((themean+qnorm(alpha)*thevariance)>1){thevariance=(1-themean)/qnorm(alpha);percentile.value=themean+qnorm(alpha)*thevariance}
    if((themean-qnorm(alpha)*thevariance)<0){thevariance=(themean)/(qnorm(alpha));percentile.value=themean+qnorm(alpha)*thevariance}
  }

  if(!is.null(themean)) {name<-"themean";value=themean;scalevalue=thevariance; percentile.value=percentile.value-0.001}
  if(!is.null(themedian)) {name<-"themedian";value=themedian; scalevalue=therange[2]-therange[1]; percentile.value=therange[2];percentile.value=percentile.value-0.001}
  if(!is.null(themode)) {name<-"themode";value=themode; scalevalue=therange[2]-therange[1]; percentile.value=therange[2];percentile.value=percentile.value-0.001}
  if (is.null(themean) && is.null(themode)){
    stopifnot ((themedian<=therange[2])|(themedian>=therange[1]))}else 
      if (is.null(themean) && is.null(themedian)){
      stopifnot ((themode<=therange[2])|(themode>=therange[1]))}
  
  
  a=runif(1,1,10)
  if (is.null(themode) && is.null(themedian))
  {to.minimize<-function(a)
  {abs(qbeta(pr_n, shape1=a, shape2=a * (1 - themean) / themean) - percentile.value)}
  }else if(is.null(themean) && is.null(themode)){
    to.minimize<-function(a)
    {abs(qbeta(pr_n, shape1=a, shape2=(3*a * (1 - themedian)+2*themedian-1) /(3*themedian)) - percentile.value)}
    }else{
      to.minimize<-function(a)
      {abs(qbeta(pr_n, shape1=a, shape2=(a * (1 - themode)+2*themode-1)/themode) - percentile.value)}}
  
  
  estimate<-optim(runif(1,1,10),to.minimize,lower = 0.001, upper = 10^4, method = "Brent")
  
  
  a_p=estimate$par
  
  if (is.null(themode) && is.null(themedian))
  { b_p = a_p * (1 - themean) / themean
  } else if (is.null(themean) && is.null(themode)){
    b_p = (3*a_p * (1 - themedian)+2*themedian-1) /(3*themedian)} else {
      b_p = (a_p * (1 - themode)+2*themode-1)/themode}
  
  set.seed(seed)
  sample_beta=rbeta(nsims,a_p,b_p)
  param=c(a=a_p,b=b_p)
  input=c(tmetric=value,scalemetric_var_or_range=scalevalue); names(input)[1]<-name
  names(input)[1]<-name
  
  if(silent==FALSE){
    cat(paste("The desired Beta distribution that satisfies the specified conditions is: Beta(", round(a_p,2), round(b_p,2),") \n",
              "Verification: The percentile value",round(qbeta(pr_n, a_p, b_p),2), "corresponds to the",pr_n,"th percentile \n"))
    cat(paste("Descriptive statistics for this distribution can be found below (or in the defined object):\n"))
    return(list(parameters=param,summary=summary(sample_beta),input=input,comment=out1))
  }
  invisible(return(list(parameters=param,summary=summary(sample_beta),input=input,comment=out1)))
  
}