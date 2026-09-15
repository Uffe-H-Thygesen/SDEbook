## Comparing the pdf with aborbing and reflecting BC's

source("fvade.R")
require(expm)

## Simple example: Brownian motion with drift
u <- function(x) -1
D <- function(x) 0.125

## Space and time
xi <- seq(0,1,length=101)
xc <- 0.5*(head(xi,-1)+tail(xi,-1))
tv <- seq(0,1,length=99)

## Generator for the absorbed and reflected case
Ga <- fvade(u,D,xi,'a')
Gr <- fvade(u,D,xi,'r')

## Start in the middle
phi0 <- numeric(length(xc))
phi0[round(length(phi0)/2)] <- 1

## Arrays for the pdfs
PHIa <- sapply(tv,function(t) phi0 %*% expm(Ga*t))
PHIr <- sapply(tv,function(t) phi0 %*% expm(Gr*t))

## Plots
require(fields)
par(mfrow=c(2,1))
image.plot(tv,xc,log(t(PHIa+1e-10)),xlab="t",ylab="x",main="P.d.f. with Asborption")
image.plot(tv,xc,log(t(PHIr+1e-10)),xlab="t",ylab="x",main="P.d.f with Reflection")
