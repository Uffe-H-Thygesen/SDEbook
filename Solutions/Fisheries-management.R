require(fields)
require(SDEtools)
require(MASS)

graphics.off()
rm(list=ls())

### Optimal fisheries management

### Parameters
K <- 1 # Carrying capacity
r <- 1 # Growth rate 
sigma <- 1/sqrt(2) # Noise in population dynamics

### Upper bound on harvest
Umax <- 16

### Discretization of state space
Xmax <- 3
dx <- 0.005
xi <- seq(0,Xmax,dx)
xc <- xi[-1] - 0.5*diff(xi)

### Functions entering into the model

### Uncontrolled system:
f <- function(x) r*x*(1-x/K)

D <- function(x) 1/2*sigma^2*x^2
dD <- function(x) sigma^2*x

u <- function(x) f(x) - dD(x)
G0 <- fvade(u,D,xi,'r')

### Effect of the fishing: The "generator" d/dx
G1 <- fvade(function(x)-1,function(x)0,xi,'r')

## Payoff
k <- function(u) c(0,sqrt(u[-1]))

uopt <- function(lambda) c(0,pmin(Umax,1/4/pmax(1e-8,-lambda[-1])^2))

sol <- PolicyIterationSingular(G0,(G1),k,(uopt),do.minimize=FALSE)

par(mfrow=c(3,1))
plot(xc,sol$V,xlab="State (biomass)",ylab="Value")
## Compare with the analytical results - note that the value function is
## only unique up to an additive constant, so shift it for better
## agreement
iref <- 10
lines(xc,0.5*log(xc) - 0.5*log(xc[iref]) + sol$V[iref])

plot(xc,sol$u,xlab="State (biomass)",ylab="Harvest")
lines(xc,xc^2)

plot(xc,sol$pi/dx,xlab="State (biomass)",ylab="P.d.f.")
lines(xc,dgamma(xc,shape=2/sigma^2-1,rate=4/sigma^2))

## We see that there is good agreement between the numerical results and the analytical ones,
## except close to the upper boundary, where the state space has been truncated.
## Fortunately, this boundary is very rarely approached, so the difference has little
## consequences.
