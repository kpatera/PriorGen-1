findbetaqq<-function(percentile.value1,percentile1,percentile.value2,percentile2)

{


findcentiles<-function(x)
{
c(
F1=qbeta(percentile1, x[1], x[2]) - percentile.value1,
F2=qbeta(percentile2, x[1], x[2]) - percentile.value2
)
}
#library(rootSolve)

ss <- multiroot(f = findcentiles, start = c(1,1))

finalshape1=ss$root[1]
finalshape2=ss$root[2]

sample_beta=rbeta(10000,finalshape1,finalshape2)

print (paste("The desired Beta distribution that satisfies the specified conditions is: Beta(", round(finalshape1,2), round(finalshape2,2),")"))
print ("Here is a plot of the specified distribution.")
#plot(density(sample_beta))
print ("Descriptive statistics for this distribution are:")
print(summary(sample_beta))
print (paste("Verification: The first percentile value",round(qbeta(percentile1, finalshape1, finalshape2),2), "corresponds to the",percentile1,"th percentile"))
print (paste("Verification: The second percentile value",round(qbeta(percentile2, finalshape1, finalshape2),2), "corresponds to the",percentile2,"th percentile"))
}
