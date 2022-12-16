install.packages("nleqslv")
require(nleqslv)
dslnex <- function(x) {
  y <- numeric(2)
  y[1] <- x[1]^2 + x[2]^2 - 2
  y[2] <- exp(x[1]-1) + x[2]^3 - 2
  y
}



xstart <- c(2,0.5)
fstart <- dslnex(xstart)
xstart
fstart
# a solution is c(1,1)
nleqslv(xstart, dslnex, control=list(btol=.01))
# Cauchy start
nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="cauchy"))
# Newton start
nleqslv(xstart, dslnex, control=list(trace=1,btol=.01,delta="newton"))
# final Broyden approximation of Jacobian (quite good)
z <- nleqslv(xstart, dslnex, jacobian=TRUE,control=list(btol=.01))
z$x
z$jac
jacdsln(z$x)


# different initial start; not a very good final approximation
xstart <- c(0.5,2)
z <- nleqslv(xstart, dslnex, jacobian=TRUE,control=list(btol=.01))
z$x
z$jac
jacdsln(z$x)
## Not run:
# no global strategy but limit stepsize
# but look carefully: a different solution is found
nleqslv(xstart, dslnex, method="Newton", global="none", control=list(trace=1,stepmax=5))
# but if the stepsize is limited even more the c(1,1) solution is found
nleqslv(xstart, dslnex, method="Newton", global="none", control=list(trace=1,stepmax=2))