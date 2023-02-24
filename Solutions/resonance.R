### Simulation of stochastic resonance

require(SDEtools)

A <- 0.12
x0 <- -1

sim <- function(sigma)
{
    Etau <- pi/sqrt(2)*exp(1/2/sigma^2)
    print(Etau)
    omega <- 2*pi/Etau
    f <- function(t,x) x-x^3 + A*cos(omega*t)
    g <- function(t,x) sigma

    times <- seq(0,min(1e5,max(10*Etau,1e3)),1e-1)
    B <- rBM(times)

    X <- heun(f,g,times,x0,B)
    X$Etau <- Etau
    X$ttau <- times/Etau 

    return(X)
}


X <- sim(0.15)

if(length(X$X)>1e4) {
    n <- length(X$X) %/% 1e5
    I <- seq(1,length(X$X),n)
    X$times <- X$times[I]
    X$ttau <- X$ttau[I]
    X$X <- X$X[I]
}

plot(X$ttau,X$X,type="l")
