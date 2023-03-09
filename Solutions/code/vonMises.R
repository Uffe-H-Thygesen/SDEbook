rm(list=ls())
graphics.off()

require(SDEtools)

xi <- seq(-pi,pi,length=101)
dx <- diff(xi)
xc <- 0.5*(tail(xi,-1)+head(xi,-1))

sigma <- 1

f <- function(x) -sin(x)
g <- function(x) sigma

u <- function(x) f(x)
D <- function(x) 0.5*sigma^2

G <- fvade(u,D,xi,'p')

rho <- StationaryDistribution(G) / diff(xi)

## Simulation

dt <- 0.01
tv <- seq(0,1000,dt)
sim <- heun(f,g,tv,0)

par(mfrow=c(1,2))
plot(tv,sim$X,type="l",xlim=c(0,50))

## Project X on the interval [-pi,pi]
Xp <- (sim$X + pi) %% (2*pi) - pi
hist(Xp,freq=FALSE)


lines(xc,rho)

Z <- 2*pi*besselI(2/sigma^2,0)
curve(1/Z*exp(sin(x)*2/sigma^2),add=TRUE,col="red")

dev.new()

acf(sin(sim$X),type="cov",lag=1000)

h <- sin(xc)
mu <- sum(rho*h*diff(xi))

acf.theory <- Vectorize(function(ti)
    sum(rho*(h-mu)*(expm(G*ti)%*% h)*diff(xi)))

tv <- seq(0,10,0.25)
lines(tv/dt,acf.theory(tv),col="red")

## Evs
require(RSpectra)

er <- eigs(G,k=4,sigma=1e-8)
el <- eigs(t(G),k=4,sigma=1e-8)

er$values <- Re(er$values)
el$values <- Re(el$values)
er$vectors <- Re(er$vectors)
el$vectors <- Re(el$vectors)

n <- nrow(G)

lambda <- sort(er$values)[3]
lines(tv/dt,acf.theory(0) * exp(lambda * tv),col="blue")

dev.new()
par(mfrow=c(3,2))
for(i in 1:3)
{
    ## Extract left eigenvector
    evl <- el$vectors[,4-i]
    ## Normalize to 2-norm = 1
    evl <- evl / sqrt(sum(evl^2/dx)) / dx

    ## Extract right eigenvector
    evr <- er$vectors[,4-i]
    evr <- evr / max(abs(evr))
    
    plot(xc,evl,type="l")
    evl2 <- evr*rho 
    evl2 <- evl2 * sum(evl^2)/sum(evl2*evl)
    lines(xc,evl2,col="red")
    plot(xc,evr,type="l",ylim=c(-1,1))
}


