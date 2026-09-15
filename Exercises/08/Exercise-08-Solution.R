require(SDEtools)
require(fields)

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
xmax <- 5
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




## Plot the stationary distribution
xmax <- 3

phi <- function(x)exp(-U(x)/D)
const <- integrate(phi,lower=-xmax,upper=xmax)
plot(function(x)phi(x)/const$value,from=-xmax,to=xmax)

tv <- seq(0,10000,0.1)
x0 <- 0

sim <- euler(f,g,tv,x0)

hist(sim$X,freq=FALSE,add=TRUE)

## Plot sample path
dev.new()
plot(sim$t[1:10000],sim$X[1:10000],type="l")


## Numerical analysis of transients
dx <- 0.01

xi <- seq(-xmax,xmax,dx)
xc <- 0.5*(head(xi,-1) + tail(xi,-1))

G <- fvade(f,function(x) D,xi,'r')

## Initial condition: A Dirac dalta at x=-1
phi0 <- 0*xc
phi0[sum(xi< -1)] <- 1

## Setup for numerical solution
T <- 100
dt <- 1

tv <- seq(dt,T,dt)

phi <- array(NA,c(length(tv),length(xc)))

## Transition probabilities
P <- expm(G*dt)

phi[1,] <- as.numeric(phi0 %*% P)

## Recursive solution of distribution
for(i in 2:length(tv))
{
    phi[i,] <- as.numeric(phi[i-1,] %*% P)
}

image.plot(tv,xc,phi)

## Backward equation
psi <- 0*t(phi)
psi[,ncol(psi)] <- as.double(xc>0)

for(i in length(tv):2)
    psi[,i-1] <- as.numeric(P %*% psi[,i])

image.plot(tv,xc,t(psi))
    
