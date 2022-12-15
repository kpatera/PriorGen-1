######################
# Exercise 1 
######################
fact <- function(){
  x <- readline(prompt = "Eisagete enan mh arnhtiko akeraio arithmo:")
  # enallaktika
  # print("Eisagete enan mh arnhtiko akeraio arithmo:") 
  # x <- scan( nmax=1)
  if(x < 0 || !(x == as.integer(x))){
    paste("O arithmos", x, "den einai mh arnhtikos akeraios", sep = " ")
  }else{
    pr <- 1
    for(i in 1:x){
      pr <- pr*i
    }
    ifelse(x == 0, paste("To paragontiko tou", x, "einai to", 1, sep = " "), 
           paste("To paragontiko tou arithmou", x, "einai to", pr))
  }
}

######################
# Exercise 2
######################
#(a) 
fx <- function( x, mu, s2, log.p=FALSE){ 
  result <- -log(x)-0.5*log(s2)-0.5*log(2*pi) - 0.5*( log(x)-mu )^2/s2
  if (!log.p) result<-exp(result)
  return(result)
}
#fx(x=1,mu=1,s2=0.01)
#fx(x=seq(0.5, 1, 0.01),mu=1,s2=0.01)


plot(fx(seq(0.5, 1, 0.01), 0.5, 1), type = 'l')
plot(seq(0.5, 1, 0.01),dlnorm(seq(0.5, 1, 0.01), 0.5, 1, log = F), type = "l") #Confirmation 
#(b)
loglike <- function( data, mu, s2 ){ 
  loglike <- sum(fx(data, mu, s2, log.p=TRUE))
  return(loglike)
}
loglike(seq(0.5, 1, 0.01), 0.5, 1)
sum(dlnorm(seq(0.5, 1, 0.01), 0.5, 1, log = T)) #Confirmation 

#(c) 
mlefx <- function(data, lowmu, himu, step=0.01){ 
  mymu <- seq(lowmu, himu, step) 
  n <- length(mymu)
  logl <- numeric(n)
  for (i in 1:n){ 
    logl[i] <- loglike( data, mu=mymu[i], s2 = 1 )
  }
  plot(logl)
  mle.mu <- mymu[which(logl == max(logl))]
  plot(mymu, logl, type='l') 
  return(mle.mu)
}		
mlefx(seq(0.5, 1, 0.01), -2, 2)