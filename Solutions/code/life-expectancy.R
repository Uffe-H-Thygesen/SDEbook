require(SDEtools)

sigma <- 0.5
mu0 <- 0.5

f <- function(x) x*(1-x)
g <- function(x) sigma*x
D <- function(x) 0.5*sigma^2*x^2
dD <- function(x) sigma^2*x
u <- function(x) f(x) - dD(x)
mu <- function(x) mu0*(1+x)

xv <- seq(0,6,0.01)
xc <- xv[-1] - 0.5*diff(xv)

G <- fvade(u,D,xv,'r') - Diagonal(x=mu(xc))
h <- solve(G,rep(-1,length(xc)))
par(mfrow=c(1,2))
plot(xc,h)
plot(xc,h*mu0)


mu0 <- 0.01
G <- fvade(u,D,xv,'r') - Diagonal(x=mu(xc))
h <- solve(G,rep(-1,length(xc)))
lines(xc,h*mu0)

mu0 <- 0.1
G <- fvade(u,D,xv,'r') - Diagonal(x=mu(xc))
h <- solve(G,rep(-1,length(xc)))
lines(xc,h*mu0)

mu0 <- 1
G <- fvade(u,D,xv,'r') - Diagonal(x=mu(xc))
h <- solve(G,rep(-1,length(xc)))
lines(xc,h*mu0)

mu0 <- 10
G <- fvade(u,D,xv,'r') - Diagonal(x=mu(xc))
h <- solve(G,rep(-1,length(xc)))
lines(xc,h*mu0)

lines(xc,1/(1+xc),col="red",lty="dashed")
