## Stochastic stability
##
## An example which is "GBM" but with additive noise so that the origin is not an equilibrium,
## and a third order term so that the system does not diverge

graphics.off()

f <- function(x) r*x - 0*abs(r)*x^3
g <- function(x) sqrt(sigma^2 *x^2 + s^2)

D <- function(x) 0.5*g(x)^2


r <- -0.5 # c(-0.6,-0.3,0,0.3,0.6)
sigmas <- c(0.6,0.7,0.9,2)
s <- 0.01

par(mfrow=c(length(sigmas),3))

for(sigma in sigmas)
{
    dx <- 0.001
    xm <- 5
    xvec <- seq(-xm,xm,dx)

    F <- function(x) f(x)/D(x)

    IF <- cumsum(F(xvec))*dx

    IF <- IF - mean(IF)

    logphiNN <- IF - log(D(xvec))

    phi <- exp(logphiNN)
    phi <- phi / sum(phi*dx)

    plot(xvec,phi,type="l")



    require(SDEtools)

    T <- 1000
    dt <- 0.01

    tv <- seq(0,T,dt)
    sim <- SDEtools::euler(f,g,tv,0)

    plot(tv,sim$X,type="l",ylim=c(-xm,xm))

    hist(pmin(pmax(sim$X,-xm),xm),breaks=xvec,main=var(sim$X))
}
