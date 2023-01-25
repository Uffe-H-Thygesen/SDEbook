## Example of the Feynmann-Kac formula being applied to compute the fitness of
## an animal

## Drift towards the surface; constant diffusivity
D <- function(x) D0
f <- function(x) -u

D0 <- 1
u <- 1

## Mortality
mu <- function(x) mu0*exp(-x/xx)
mu0 <- 1
xx <- 2

## Feeding
h <- function(x) 1/(1+exp(x-2))
    
## maximum Depth computed from the diffusive boundary layer
H <- 5*D0/u
xi <- seq(0,H,length=101)
xc <- 0.5*(head(xi,-1)+tail(xi,-1))
    
require(SDEtools)

## Compute generator
G <- fvade(f,D,xi,'r')

## Subtract mortality from the generator to get the Feynmann-Kac generator 
G <- G - Diagonal(length(xc),x=mu(xc))

## Compute fitness
V <- solve(G,-h(xc))

plot(V,xc,ylim=c(H,0),type="l")
