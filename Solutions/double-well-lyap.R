## Compute the Lyapunov exponent for the double well system

require(SDEtools)

## Here we show the results only for this parameter value:
mu <- 1

## Dynamics
f <- function(x) mu*x - x^3
g <- function(x) 0.5

## Jacobians
fp <- function(x) mu - 3*x^2
gp <- function(x) 0

## Combined dynamics of state and sensitivity
ff <- function(xs) c(f(xs[1]),fp(xs[1]) * xs[2])
gg <- function(xs) c(g(xs[1]),gp(xs[1]) * xs[2])

T <- 200
dt <- 0.1

times <- seq(0,T,dt)

xs0 <- c(0,1)

sim <- heun(ff,gg,times,xs0)

par(mfrow=c(2,1))
plot(times,sim$X[,1],type="l",xlab="Time",ylab="X")
plot(times,log(sim$X[,2])/times,type="l",xlab="Time",ylab=expression(mu[t]))

tail(log(sim$X[,2])/times)

## Potential, its double derivative, and the unnormlized p.d.f.
D <- 0.5*sigma^2
U <- function(x) -0.5*mu*x^2  + 0.25*x^4
Upp <- function(x) -mu + 3*x^2
phi <- function(x) exp(-U(x)/D)

Z <- integrate(phi,lower=-Inf,upper=Inf)$value
lambda <- -integrate(function(x)phi(x)*Upp(x),lower=-Inf,upper=Inf)$value / Z

abline(h=lambda,lty=2)
legend("topright",legend=c("Simulation","Numerical"),lty=1:2)
