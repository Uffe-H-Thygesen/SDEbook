## Eigenvalues of the Kolmogorov equations for stochastic logistic growth

require(SDEtools)
require(RSpectra)

f <- function(x) x*(1-x)
g <- function(x) s*x

D <- function(x) s^2*x^2/2
Dp <- function(x) s^2 * x
u <- function(x) f(x) - Dp(x)

xi <- seq(0,4,0.001)
xc <- 0.5*(head(xi,-1)+tail(xi,-1))

s <- 1.4

G <- fvade(u,D,xi,'r')
n <- nrow(G)

# evs <- eigen(G)

evs <- eigs(G,2,"LM",sigma=0.01)


print(max(-Re(evs$values)))
## plot(xc,Re(evs$vectors[,1]))

## Lyap
-(1-s^2/2)
