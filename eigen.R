## Example of eigenfunctions: The Rosenzweig-MacArthur example

require(SDEtools)

fII <- function(N) beta*N/(beta*N/Cmax + 1)
fN <- function(NP) r*NP[1]*(1-NP[1]/K) - fII(NP[1])*NP[2]
fP <- function(NP) epsilon*fII(NP[1])*NP[2] - mu * NP[2]

f <- function(NP) c(fN(NP),fP(NP))
g <- function(NP) diag(NP*sigmaNP)

## Parameters
K <- 1
r <- 1
sigmaNP <- c(0.05,0.05)
beta <- 5
Cmax <- 2
epsilon <- 0.125
mu <- 0.1

NP0 <- c(0.1,0.01)

times <- seq(0,1000,0.01)
sol <- euler(f,g,times,NP0)

matplot(sol$times,sol$X,type="l")


nv <- seq(0,2,0.04)
pv <- seq(0,1,0.025)

nc <- 0.5*(head(nv,-1)+tail(nv,-1))
pc <- 0.5*(head(pv,-1)+tail(pv,-1))

nn <- length(nc)
np <- length(pc)

D <- function(NP) 0.5*g(NP)^2
u <- function(NP) f(NP) - sigmaNP^2 * NP

G <- fvade2d(function(N,P)u(c(N,P))[1],
             function(N,P)u(c(N,P))[2],
             function(N,P)D(c(N,P))[1,1],
             function(N,P)D(c(N,P))[2,2],
             nv,pv)

rho <- StationaryDistribution(G)

rho <- unpack.field(rho,nn,np)
image(nc,pc,t(rho))

points(sol$X[,1],sol$X[,2],pch=".")
