## Black Scholes formula - comparing the analytical result vs a numerical solution

require(SDEtools)

## Parameters
r <- 0.05
sigma <- 0.02
K <- 1

dt <- 0.1
T <- 10 
tvec <- seq(0,T,dt)

## Numerical grid. This grid extends longer than what we aim
## to plot, so that numerical boundary effects do not appear
## in the plot
xvec <- seq(0,5,0.01)
xc <- 0.5*(head(xvec,-1)+tail(xvec,-1))

## The analytical result
price <- function(x,t)
{
    d <- (log(x/K) + (r-0.5*sigma^2)*(T-t))/sigma/sqrt(T-t)
    x*pnorm(d+sigma*sqrt(T-t)) - exp(-r*(T-t))*K*pnorm(d)
}

## Tabulate the price
V <- outer(xvec,tvec,price)

## Numerical solution. 
## Advection-diffusion formalism
u <- function(x) r*x-sigma^2*x
D <- function(x) 0.5*sigma^2*x^2                     
G <- fvade(u,D,xvec,'r')
eGt <- expm(G*dt)

V1 <- array(NA,c(length(xc),length(tvec)))
V1[,length(tvec)] <- pmax(0,xc-K)
ert <- exp(r*dt)
for(i in length(tvec):2) V1[,i-1] <- as.numeric(eGt %*% V1[,i]) / ert

par(mfrow=c(1,2))
n <- sum(xvec < 2)
levels <- c(10^seq(-7,-1,1),seq(0.2,1.3,0.1))
contour(tvec,xvec[1:n],t(V[1:n,]),levels=levels,main="Analytical")
contour(tvec,xc[1:n],t(V1[1:n,]),levels=levels,main="Numerical")

