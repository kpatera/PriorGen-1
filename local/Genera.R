require(devtools)
require(remotes)
remotes::install_github("kpatera/PriorGen-1",force = T)
library(PriorGen)

out<-findbeta(themode=0.15, percentile=0.90,lower.v=TRUE, percentile.value=0.40)
out2<-findbeta(themode=0.20, percentile=0.90,lower.v=TRUE, percentile.value=0.40)
plot(out,type="l")
print(out)

plot(out, type="l",lwd=3)
plot(out, type="l",lwd=3,col="red")

temp=findbetamupsi(themean=0.10, percentile=0.79, 
                   lower.v=TRUE, percentile.value=0.26, 
                   psi.percentile=0.95, percentile.median=0.28, 
                   percentile95value=0.3)
temp$param_gamma
temp$param_beta

N(logit(0.2), 0.3 2)
LO c + 𝜃 + 𝜏Z

plot(density(rnorm(n = 10000,mean = log(0.2/0.8),sd = 0.3^2)),xlim=c(-2,2))
lines(density(rnorm(n = 10000,mean = log(0.2/0.8)+1+0+0.05*rnorm(10000),sd = 0.3^2)))
