findbetaqq<-function(percentile.value1,percentile1,percentile.value2,percentile2,
                     seed=280385, silent=TRUE){
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
  
  if(silent==FALSE){
    print (paste("The desired Beta distribution that satisfies the specified conditions is: Beta(", round(finalshape1,2),",", round(finalshape2,2),")"))
    print ("Descriptive statistics for this distribution are:")
    #print(summary(sample_beta))
    (paste("Verification: The first percentile value",round(qbeta(percentile1, finalshape1, finalshape2),2), "corresponds to the",percentile1,"th percentile"))
    (paste("Verification: The second percentile value",round(qbeta(percentile2, finalshape1, finalshape2),2), "corresponds to the",percentile2,"th percentile"))
    # print(paste("The procedure took",round((proc.time()-time)[3],3),"minutes and ",k ,"loops"))
    return(list(parameters=param,summary=summary(sample_beta),input=input))
  }
  invisible(return(list(parameters=param,summary=summary(sample_beta),input=input)))
}
