require(SDEtools)

## CIR
f <- function(x) lambda*(xi-x)
g <- function(x) gamma*sqrt(x)

## Lamperti transformed, Z=sqrt(X)
fL <- function(z) f(z^2)/2/z - 0.1255*g(z^2)/z^3
gL <- function(z) g(z^2)/2/z

lambda <- 1
xi <- 1
gamma <- 1

T <- 10
dt <- 0.01

tv <- seq(0,T,dt)

B <- rBM(tv)

x0 <- xi

solE <- euler(f,g,tv,x0=c(X=x0),B=B,p=abs)

plot(solE$t,solE$X,type="l")

solL <- euler(fL,gL,tv,x0=c(Z=sqrt(x0)),B=B,p=abs)

solL$Z <- solL$X
solL$X <- solL$Z^2


lines(solL$t,solL$X,col="red")
