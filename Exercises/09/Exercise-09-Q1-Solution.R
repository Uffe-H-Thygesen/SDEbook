## Comparing the pdf with aborbing and reflecting BC's

require(SDEtools)

## Simple example: Brownian motion with drift
u <- function(x) -1
D <- function(x) 0.125

## Space and time
xi <- seq(0,1,length=101)
xc <- 0.5*(head(xi,-1)+tail(xi,-1))
dx <- diff(xi)
tv <- seq(0,5,length=99)

## Generator for the absorbed and reflected case
Ga <- fvade(u,D,xi,'a')
Gr <- fvade(u,D,xi,'r')

## Start in the middle
phi0 <- numeric(length(xc))
phi0[round(length(phi0)/2)] <- 1

## Arrays for the pdfs
PHIa <- sapply(tv,function(t) as.numeric(phi0 %*% expm(Ga*t)))
PHIr <- sapply(tv,function(t) as.numeric(phi0 %*% expm(Gr*t)))

## Plots
require(fields)
par(mfrow=c(2,1))
image.plot(tv,xc,log(t(PHIa+1e-10)),xlab="t",ylab="x",main="P.d.f. with Asborption")
image.plot(tv,xc,log(t(PHIr+1e-10)),xlab="t",ylab="x",main="P.d.f with Reflection")


## Compute eigenvalue with largest real part
evs <- eigen(t(Ga))  # Note: We could have used a sparse library (spam RSpectra)
i <- which.max(Re(evs$values))

lambda <- evs$values[i]
phi <- evs$vectors[,i]
phi <- phi/sum(phi*dx)

plot(xc,phi,type="l")
points(xc,PHIa[,length(tv)] / sum(PHIa[,length(tv)] * dx),col="blue")

pi <- StationaryDistribution(Gr)
lines(xc,pi/sum(pi*dx),lty="dashed")
