## Application of the Feynmann-Kac formula to compute the present-day value of a future
## under the CIR model for the interest rate

## CIR parameters
xi <- 0.05  ## Mean interest rate (per year)
lambda <- 1 ## Relaxation rate of interest rate
gamma <- 0.2 ## Volatility

## Drift, diffusion, advection
f <- function(x) lambda*(xi-x)
D <- function(x) 0.5*gamma^2*x
dD <- function(x) 0.5*gamma^2 + 0*x
u <- function(x) f(x) - dD(x)

require(SDEtools)
xv <- seq(0,4*xi,length=101)
xc <- 0.5*(head(xv,-1)+tail(xv,-1))

G <- fvade(u,D,xv,'r')

pi <- StationaryDistribution(G)

par(mfrow=c(2,1))
plot(xc,pi/diff(xv))

M <- G - Diagonal(length(xc),x=xc)

T <- 10
psi <- rowSums(expm(M*T))
plot(xc,psi)
