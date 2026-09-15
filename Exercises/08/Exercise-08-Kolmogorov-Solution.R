## ----message=FALSE------------------------------------------------------------
require(SDEtools)
require(fields)


## -----------------------------------------------------------------------------
## Define dynamics

## Drift and intensity
f <- function(x) lambda * (xi - x)
g <- function(x) gamma*sqrt(abs(x))

## Diffusivity and its spatial derivative
D <- function(x) 0.5*gamma^2*x
Dp <- function(x) 0.5*gamma^2

## Advective flow field
u <- function(x) f(x) - Dp(x)

## Parameters
xi <- 2
gamma <- 1
lambda <- 1/2

## Spatial grid
xmax <- 10
xv <- seq(0,xmax,length=101)
dx <- diff(xv)
xc <- xv[-1] - 0.5*dx

## Discretize the generator
G <- fvade(u,D,xv,'r')

## Compute the stationary density
pi <- StationaryDistribution(G)
phi <- pi / dx
plot(xc,phi,type="l")
plot(function(x)dgamma(x,rate=2*lambda/gamma^2,shape=2*lambda*xi/gamma^2),
     from=0,to=xmax,add=TRUE,col="red")


## -----------------------------------------------------------------------------
## Transients

## Initial condition for the SDE
x0 <- xi/4

## Initial condition for the FKE
phi0 <- numeric(length(xc))
phi0[sum(xc<x0)] <- 1

## Time grid
tv <- seq(0,10,0.1)

## Solve the FKE
PHI <- sapply(tv,function(t) as.numeric(phi0 %*% expm(G*t)))/dx

image.plot(tv,xc,t(PHI))

CDF <- apply(PHI*dx,2,cumsum)

image.plot(tv,xc,t(CDF))


## -----------------------------------------------------------------------------
plot(xc,PHI[,length(tv)])
lines(xc,phi)


## -----------------------------------------------------------------------------
  EX <- apply(PHI*dx*xc,2,sum)

  EX2 <- apply(PHI*dx*xc^2,2,sum)
  VX <- EX2 - EX^2
  sX <- sqrt(VX)

  plot(tv,EX,lwd=2,ylim=c(0,max(EX+sX)))
  lines(tv,EX+sX,lty="dashed")
  lines(tv,EX-sX,lty="dashed")

  EXa <- xi+(x0-xi)*exp(-lambda*tv)
  lines(tv,EXa,col=2)

  sXa <- sqrt(gamma^2*xi/2/lambda)

  abline(h=xi,col=2,lty=2)
  abline(h=xi+sXa,col=2,lty=2)
  abline(h=xi-sXa,col=2,lty=2)


## -----------------------------------------------------------------------------
  h <- (xc>=2)
  T <- tail(tv,1)
  k <- sapply(tv,function(t) as.numeric(expm(G*(T-t)) %*% h))
  image.plot(tv,xc,t(k))


## -----------------------------------------------------------------------------
tsim <- seq(0,1000,0.1)

sim <- euler(f,g,tsim,x0,p=abs)
hist(sim$X,freq=FALSE,breaks=50)
lines(xc,phi)

