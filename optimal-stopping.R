## Example of optimal stopping


f <- function(x) 0
g <- function(x) 1

s <- function(x) sin(x)

xi <- seq(0,2*pi,length=101)
xc <- 0.5*(head(xi,-1)+tail(xi,-1))

sc <- s(xc)

require(SDEtools)

G <- fvade(f,g,xi,'r')

uopt <- function(Wp) as.numeric(Wp>0)

k <- function(u) G1 %*% u[[1]]

## Policyiteration

## Start with an arbitrary strategy
u <- as.numeric(runif(length(xc))>0.5)

plot(xc,sc,type="l")

for(i in 1:2)
{
    G0 <- diag(1-u)
    G1 <- G * as.numeric(u)
    V <- as.numeric(solve(G0+G1,G0%*%sc))
    u <- as.numeric((G%*%V) >= 0 )
    lines(xc,V)
}

lines(xc,V,lwd=3)
