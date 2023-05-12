#' The findbetamodel function (UNDER CONSTRUCTION)
#'
#' A function to estimate the parameters of a Beta distribution for the expected mean of a proportion.
#' This function first produces a posterior MCMC-based distribution for the prevalence based on a meta-analysis of the available studies.
#' Then, a Beta(a,b) distribution is elicited based on the available the produced MCMC samples.
#'
#' @usage findbetamutau(pos.c, n.c,
#'   seed = 280385, nsims = 10000, root.method = "multiroot")
#'
#' @param pos.c specify a vector of positive cases for each study. It takes a values between 0 and 1.
#' @param n.c specify a vector of the total sample size for each study. It takes a values between 0 and 1.
#' @param seed A fixed seed for replication purposes.
#' @param nsims Number of simulations for the creation of various summary metrics of the elicited prior.
#' @param root.method Choose between two alternatives to solve the two non-linear equations to identify the hyperparameters of psi. root.method="multiroot" involves the basic function of the rootSolve package, root.method="nleqslv" involves the base functions of the nleqslv package.
#'
#' @examples
#' ## Example
#' ## The mean prevalence of a disease/infection for the units within an area/region
#' ## is thought to be 0.20 and we are 99% confident that it is not more than 0.40.
#' ## Within this area, we are also confident that 90% of all units have a prevalence
#' ## less or equal to 0.50 and we are 95% certain that it does not exceed 0.60
#'
#' res_meta <- findbetamodel(
#'   pos.c = c(20,10,30,13,15,8,24), n.c = c(1000,1000,900,800,1100,600,900)
#' )
#'
#' @export
#' @return param_beta: The beta distribution parameters Beta(a,b)
#' @return param_gamma: The gamma distribution parameters Gamma(a,b)
#' @return summary: A basic summary of the elicited prior
#' @return input: NULL
#' @return param_upper: NULL
#'
#' @references
#' 
#' 
findbetamodel <- function(pos.c, n.c, seed = 280385, 
                          nsims = 10000, root.method = "multiroot") {
  require(rjags)
  require(MASS)
  cat("model{

for(i in 1:k){
y[i] ~ dbin(ap[i], n[i])
ap[i] <- sub.p[i]*main.Se + (1-sub.p[i])*(1-main.Sp)

sub.p[i] ~ dbeta(alpha,beta) T(0.001,0.999)
}

alpha <- main.ap*main.psi
beta <- main.psi*(1-main.ap)
main.ap ~ dbeta(aa, bb) T(0.001,0.999)
main.psi ~ dgamma(0.01, 0.01) T(0.001,0.999)

#informative prior for Se and Sp
main.Se <- 1
main.Sp <- 1

#predictions
y.pre ~ dbin(main.pstar.rep,m)

main.pstar.rep ~ dbeta(alpha,beta)
pre.pequal0 <- equals(main.pstar.rep,0)
pre.plessthan0.05 <- step(0.05-main.pstar.rep)
}", file=paste("findbetamodel.txt"))
  
  
  
  SaveParams <- c("main.ap","main.psi","main.pstar.rep","sub.p","pre.pequal0",
                  "pre.plessthan0.05","y.pre") 
  
  
  
  
  jagsoutput_AppMult<-rjags::jags.model(data=list(n=n.c,y=pos.c,m=100,
                                                  k=length(n.c),
                                                  aa=1, bb=1),
                                        inits=NULL, n.chains=1,
                                        n.adapt = floor(nsims/10),
                                        file=paste("findbetamodel.txt"),quiet=TRUE)
  
  model1 <<- coda.samples(jagsoutput_AppMult,
                               n.iter=nsims,thin = 3,
                               variable.names=SaveParams,seed=998)
  
  param_beta<-fitdistr(x = as.numeric(model1[[1]][,1]),densfun = "beta",start = list(shape1=1,shape2=1), lower=c(0,0))
  param_gamma<-fitdistr(x = as.numeric(model1[[1]][,2]),densfun = "gamma",start = list(shape=1,rate=1), lower=c(0,0))
  
  out <- list(param_beta = param_beta, param_gamma = param_gamma, summary = summary(model1[[1]][,3]), 
              input = NA, param_upper = NA)
  class(out) <- "PriorGen"
  return(out)
}
