### Sample path of Ornstein-Uhlenbeck process.
### Sample paths as well as 1-sigma limits

### Note: In stead of simulating the process with the Euler-Maruyama
### method (or similar), we simulate the Ito integrals and transform,
### following the analytical derivation of the solution formula.

set.seed(123456)
require(SDEtools)

## Parameters in the model dX=A*X*dt + G*dB
A <- -1
G <- 1

## Simulation control
T <- 2
h <- 0.001
tvec <- seq(0,T,h)

B <- rBM(tvec)

## Solve for Y=exp(-A*t)*X, which satisfies dY = G*exp(-A*t)*dB 
Y <- itointegral(G*exp(-A*tvec),B)

## transform to X=exp(A*t)*Y
Xt <- exp(A*tvec)*Y

## Find the variance of Y from the Ito isometry
VY <- itointegral(G^2*exp(-2*A*tvec),tvec)

## Transform to the variance of X
VXt <- exp(2*A*tvec) * VY

pdf("ito-integral-linear.pdf")

plot(tvec,Xt,type="l",ylim=c(-1,1)*2*sqrt(VXt[length(tvec)]),
     xlab="t",ylab=expression(X[t]))
lines(tvec,sqrt(VXt),lty="dashed")
lines(tvec,-sqrt(VXt),lty="dashed")

dev.off()
