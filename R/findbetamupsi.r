findbetamupsi<-function(themean, percentile=0.95, lower.v=T, percentile.value, 
                        psi.percentile=0.90, percentile.median, percentile95value,
                        seed=280385,silent=TRUE){
  if (lower.v == T){pr_n = percentile}else{pr_n = 1 - percentile}
  stopifnot((lower.v == T && themean <= percentile.value) | 
              (lower.v == F && themean >= percentile.value))
  stopifnot(percentile.median > themean && percentile95value > 
              percentile.median)
  stopifnot((lower.v == T && percentile.value < percentile.median) | 
              (lower.v == F && percentile.value != percentile.median))
  stopifnot(psi.percentile > 0.5 && psi.percentile < 1)
  a = runif(1, 1, 10)
  to.minimize <- function(a) {
    abs(qbeta(pr_n, shape1 = a, shape2 = a * (1 - themean)/themean) - 
          percentile.value)
  }
  estimate <- optim(runif(1, 1, 10), to.minimize, lower = 0.1, 
                    upper = 10^4, method = "Brent")
  alpha_mu = estimate$par
  beta_mu = alpha_mu * (1 - themean)/themean
  mu = alpha_mu/(alpha_mu + beta_mu)
  gu = alpha_mu + beta_mu
  f <- function(ga) abs(qbeta(psi.percentile, mu * ga, ga * 
                                (1 - mu)) - percentile.median)
  r <- optim(1, f, lower = 0, upper = 10^4, method = "Brent")
  alpha_t = r$par[1]
  f <- function(ga) abs(qbeta(psi.percentile, mu * ga, ga * 
                                (1 - mu)) - percentile95value)
  r <- optim(1, f, lower = 0, upper = 10^4, method = "Brent")
  beta_t = r$par[1]
  model <- function(x) c(F1 = qgamma(0.5, shape = x[1], scale = 1/x[2]) - 
                           alpha_t, F2 = qgamma(0.05, shape = x[1], scale = 1/x[2]) - 
                           beta_t)
  ss2 <- multiroot(f = model, start = c(5, 2), positive = T)
  qgamma(0.5, ss2$root[1], ss2$root[2])
  qgamma(0.95, ss2$root[1], ss2$root[2])
  
  
  a=rgamma(10000,ss2$root[1],ss2$root[2])
  b=rbeta(10000,alpha_mu,beta_mu)
  param=c(a=alpha_mu,b=beta_mu)
  sample_beta=rbeta(10000,a*b,a*(1-b))
  
  input=c(themean=themean, percentile=percentile,
          percentile.value=percentile.value, psi.percentile=psi.percentile, 
          percentile.median=percentile.median, percentile95value=percentile95value)
  
  if(silent==FALSE){
    print (paste("The desired Beta distribution that satisfies the specified conditions about the mean of the prevalence 'mu' is: Beta(", round(finalshape1,2), round(finalshape2,2),")"))
    print (paste("The desired Gamma distribution that satisfies the specified conditions about the variability 'psi' of the prevalence is: Gamma(", round(ss2$root[1],2), round(ss2$root[2],2),")"))
    #print ("The plot gives the specified prior beleif on the prevalence distribution.")
    #plot(density(rbeta(10000,a*b,a*(1-b))))
    print("Descriptive statistics for this distrubiton are:")
    return(list(parameters=param,bot_param=list(at=a*b,bt=a*(1-b)),summary=summary(sample_beta),input=input))
  }
  invisible(return(list(parameters=param,bot_param=list(at=a*b,bt=a*(1-b)),summary=summary(sample_beta),input=input)))
  
}
