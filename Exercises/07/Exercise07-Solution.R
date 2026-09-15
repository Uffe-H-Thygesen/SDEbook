## Start from a clean slate
rm(list=ls())
graphics.off()
require(SDEtools)

## Define drift and noise terms
f <- function(x) c(x[2],sin(x[1]) - lambda*x[2])
g <- function(x) c(0,sigma)

## System parameters
lambda <- 0.1
sigma <- 0.1

## Simulation parameters
dt <- 0.01
x0 <- c(pi,sqrt(0.5)*sigma)
T <- 10000

tv <- seq(0,T,dt)

## Perform the simulation with the Euler method
B <- rBM(tv)
sol <- euler(f,g,tv,x0,B)

par(mfrow=c(1,2))
plot(sol$times,sol$X[,1],type="l")
hist(sol$X[,1])

## The Ito equation for the energy is dE = (-lambda*V^2 + 1/2*sigma^2)*dt  + V*sigma*dB
## In steady-state, the drift term must have expectation zero
## So E(V^2/2) = sigma^2/(4*lambda)

BurnIn <- round(0.1*length(tv)):length(tv)

Ek <- 0.5*sol$X[,2]^2
print(c(EkEmp=mean(Ek[BurnIn]),EkAnal=sigma^2/(4*lambda)))

sol2 <- heun(f,g,tv,x0,B)

Ek2 <- 0.5*sol2$X[,2]^2
print(c(EkEmp=mean(Ek2[BurnIn]),EkAnal=sigma^2/(4*lambda)))

dev.new()
I <- 1:10000
plot(sol$times[I],sol$X[I,1],type="l")
lines(sol2$times[I],sol2$X[I,1],col="red")
