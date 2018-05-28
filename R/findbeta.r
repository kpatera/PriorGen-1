

findbeta<-function(themean=NULL,themedian=NULL,themode=NULL,percentile=0.95, lower.v=F, percentile.value)

{
stopifnot ( (is.null(themean) + is.null(themedian) + is.null(themode)) == 2 )


  if (is.null(themode) && is.null(themedian))
  {stopifnot ((lower.v==T && themean<=percentile.value)|(lower.v==F && themean>=percentile.value))
  }else if (is.null(themean) && is.null(themode)){
    stopifnot ((lower.v==T && themedian<=percentile.value)|(lower.v==F && themedian>=percentile.value))}else{
      stopifnot ((lower.v==T && themode<=percentile.value)|(lower.v==F && themode>=percentile.value))}


a=runif(1,1,10)

if (lower.v==T){pr_n=percentile} else {
pr_n=1-percentile}




if (is.null(themode) && is.null(themedian))
{to.minimize<-function(a)
{abs(qbeta(pr_n, shape1=a, shape2=a * (1 - themean) / themean) - percentile.value)}
}else if (is.null(themean) && is.null(themode)){
to.minimize<-function(a)
{abs(qbeta(pr_n, shape1=a, shape2=(3*a * (1 - themedian)+2*themedian-1) /(3*themedian)) - percentile.value)}}else{
to.minimize<-function(a)
{abs(qbeta(pr_n, shape1=a, shape2=(a * (1 - themode)+2*themode-1)/themode) - percentile.value)}}


estimate<-optim(runif(1,1,10),to.minimize,lower = 0.1, upper = 10^4, method = "Brent")


finalshape1=estimate$par

if (is.null(themode) && is.null(themedian))
{ finalshape2 =finalshape1 * (1 - themean) / themean
} else if (is.null(themean) && is.null(themode)){
finalshape2 = (3*finalshape1 * (1 - themedian)+2*themedian-1) /(3*themedian)} else {
finalshape2= (finalshape1 * (1 - themode)+2*themode-1)/themode}

print (paste("The desired Beta distribution that satisfies the specified conditions is: Beta(", round(finalshape1,2), round(finalshape2,2),")"))
sample_beta=rbeta(10000,finalshape1,finalshape2)
print ("Here is a plot of the specified distribution.")
#plot(density(sample_beta))
print ("Descriptive statistics for this distribution are:")
print(summary(sample_beta))
print (paste("Verification: The percentile value",round(qbeta(pr_n, finalshape1, finalshape2),2), "corresponds to the",pr_n,"th percentile"))

}


