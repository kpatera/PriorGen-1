require(devtools)
require(remotes)
remotes::install_github("kpatera/PriorGen-1",force = T)
#library(PriorGen)

out<-findbeta(themode=0.15, percentile=0.90,lower.v=TRUE, percentile.value=0.40)
plot(out,type="l")
print(out)

findbeta_plot(out, type="l",lwd=3,)
