# For figure illustrating the reflection argument for the maximum of Brownian motion

pdf(file="BMreflection.pdf")

set.seed(12)  # For reproducible results

## Time grid
dt <- 0.001
T <- 10
tv <- seq(0,T,dt)
Nt <- length(tv)

## Simulate Brownian motion 
dB <- sqrt(dt)*rnorm(length(tv)-1)
B <- c(0,cumsum(dB))

## Threshold
x <- 2

## Identify the first index beyond x 
tauI <- which(B>=x)[1]

## Set up reflected BM
Br <- B

## If hit, then reflect
if(is.finite(tauI)) Br[(tauI+1):Nt] <- 2*x - B[(tauI+1):Nt]

## Plotting ranges
Blim <- range(c(B,Br))

par(cex=1.5)
plot(tv,Br,type="l",lwd=3,ylim=Blim,xlab="t",ylab=expression(B[t]))
lines(tv,B,lty="dotted",lwd=1)
lines(c(0,T),rep(x,2),lty="dotted")

dev.off()
