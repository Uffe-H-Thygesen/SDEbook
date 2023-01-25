### Figure comparing Euler vs. the analytical solution of of Geometric Brownian Motion

set.seed(1111111)

rm(list=ls())
graphics.off()

## GBM parameters
r <- 0.5
sigma <- 0.5
x <- 1

## Simulation parameters
T <- 1
dt <- 1e-4 # For plot of analytical solution
h <- 0.25  #T ime step for Euler (should divide dt)

tfine <- seq(0,T,dt)
tcoarse <- seq(0,T,h)

require(SDEtools)

Bfine <- rBM(tfine)
Bcoarse <- approx(tfine,Bfine,tcoarse)$y 

## Analytical solution
Xa <- x*exp((r-0.5*sigma^2)*tfine + sigma*Bfine)

## Euler-Maruyama  approximation (analytical expression)
Xe <- x*c(1,cumprod(1+r*h+sigma*diff(Bcoarse)))

pdf(file="figEuler.pdf",width=5,height=4)

plot(tfine,Xa,type="l",ylim=range(c(Xe,Xa)),col="grey",xlab="Time",ylab="Solution X")
lines(tcoarse,Xe)
legend("topleft",c("Analytical","Euler"),col=c("grey","black"),lty="solid",pch=c(NA,1))

points(tcoarse,Xe,pch=1)

dev.off()
