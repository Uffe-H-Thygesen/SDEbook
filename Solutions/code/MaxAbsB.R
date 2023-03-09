### Find the distribution of the maximum of Brownian motion
###
### uht, August 9, 2011

##########################################################
## Monte Carlo simulation
N <- 500              # Number of sample paths
T <- 1                # Final time
dt <- 0.001           # Time step

nt <- round(T/dt)     # Number of time steps
tvec <- seq(0,T,dt)   # Vector of time points. Note: Length nt+1

## Simulate Brownian motion. 
dB <- array(sqrt(dt)*rnorm(N*nt),c(nt,N))
B <- apply(dB,2,function(x)c(0,cumsum(x)))

## Helper function to compute the running max of a path
cummaxabs <- function(x)
  {
    cma <- numeric(length(x))
    cma[1] <- abs(x[1])
    for(i in 2:length(x)) cma[i] <- max(abs(x[i]),cma[i-1])
    return(cma)
  }

## Compute running max for each path
S <- apply(B,2,cummaxabs)

graphics.off()

## Compute quantiles of the S process, at each time step
quants <- t(apply(S,1,quantile))

## Set up plot of S
plot(c(0,T),c(0,max(quants[nt,4],S[nt,1])),
              type="n",xlab="t",
              ylab=expression(S[t]))

## Grey shaded region to include confidence regions
polygon(c(tvec,tvec[(nt+1):1]),
        c(quants[,2],quants[(nt+1):1,4]),
        border=NA,
        col="grey")

## Solid thin line to show one sample path
lines(tvec,S[,1],lwd=1)

## Thick solid line to show the median
lines(tvec,quants[,3],lwd=3)

dev.copy2pdf(file="maxAbsB-sim.pdf")

## Setup new plot: The distribution of S[1]
X11()
plot(sort(S[nt+1,]),(1:N)/N,type="l",lwd=2,
     xlab="s",
     ylab=expression("Empirical c.d.f. of "*S[1]*", "*F[S[1]](s)))

## Save plot for later, to ammend
DistPlot <- dev.cur()

X11()

## Tabulate "survival function", i.e. P(no absorption yet)
SurvTab <- apply(S<1,1,mean)

plot(tvec,SurvTab,type="l",
     xlab="t",ylim=c(0,1),
     lwd=1,
     ylab=expression(P(tau>t)))

SurvPlot <- dev.cur()

## Add estimate of c.d.f., derived from the survival function
## using the scaling properties of S
dev.set(DistPlot)
lines(1/sqrt(tvec),SurvTab,col="red",lwd=2)

##################################################
## Finite volume discretization

# Grid of [-1,1]: Number of grid cells and cell width
nx <- 101
dx <- 2/nx

## Vector containing centers of grid cells, for plotting
xc <- seq(-1+dx/2,1-dx/2,length=nx)

## The number of the cell containing the origin
n0 <- ceiling(nx/2) 
  
## Construct discretized generator
G <- array(0,c(nx,nx))
diag(G) <- rep(-1/dx^2,nx)
for(i in 2:nx)
  {
    G[i,i-1] <- 0.5/dx^2
    G[i-1,i] <- 0.5/dx^2
  }

## Load function to compute matrix exponential
source("expm.R")

## Tabulate "survival function", i.e. probability
## that the process has stayed within (-1,1) so far

## Setup time points and empty array for survival function
ts <- seq(0,sqrt(10),0.1)^2  # increased resolution at small t
Gt <- numeric(length(ts))

for(i in 1:length(ts))
  {
    ## Compute transition probability
    Pt <- expm(G*ts[i])

    ## Distribution when starting at the origin
    phi <- Pt[n0,]

    ## Sum to find survival probability
    Gt[i] <- sum(phi)
  }

## Ammend graph to empirical survival function
dev.set(SurvPlot)
lines(ts,Gt,lwd=2,col="blue")

dev.copy2pdf(file="Absorbed-BM-survival.pdf")

## ... rescale, and add to c.d.f. plot
dev.set(DistPlot)
lines(1/sqrt(ts),Gt,col="blue",lwd=2)

dev.copy2pdf(file="Time-to-absorption-cdf.pdf")

### Set up plot of density function at specified times

tplot <- c(0.05,0.1,1)

Pt <- expm(G*tplot[1])
phi <- Pt[n0,]/dx

X11()
plot(xc,phi,col="blue",type="l",
     lty=1,
     lwd=4,xlab="x",
     ylab=expression(phi(x,t)))

DensPlot <- dev.cur()

for(i in 2:length(tplot))
  {
    Pt <- expm(G*tplot[i])
    phi <- Pt[n0,]/dx
    lines(xc,phi,lwd=4,col="blue",lty=i)
  }

### Compute densities with Fourier method

xvec <- seq(-1,1,0.01) ## Grid, for plotting
nf <- 10               ## No. Fourier coefficients

### Helper function to compute the solution at time t using truncated
### Fourier series.

FourierSolution <- function(t)
  {
    ## Initialize phi
    phiF <- numeric(length(xvec))

    ## Add each term in the Fourier series
    for(i in 0:nf)
      {
        k <- 2*i +1
        phiF <- phiF + exp(-pi^2/8*k^2*t)*cos(pi/2*k*xvec)
      }

    return(phiF)
  }

for(i in 1:length(tplot))
  {
    phiF <- FourierSolution(tplot[i])
    
    ## Add to plot
    lines(xvec,phiF,lwd=1,lty=1,col="red")
  }

dev.copy2pdf(file="Absorbed-BM-density.pdf")


### Plot of periodic extension
X11()

phiF <- FourierSolution(0.1)
plot(xvec,phiF,type="l",lwd=4,xlab="x",
     xlim=c(-4,4),ylim=c(-1,1)*max(phiF),
     ylab=expression(phi(x)))

# Extensions
lines(xvec-4,+phiF,lwd=3,col="red")
lines(xvec-2,-phiF,lwd=3,col="red")
lines(xvec+2,-phiF,lwd=3,col="red")
lines(xvec+4,+phiF,lwd=3,col="red")

lines(xvec-2,+phiF,lwd=3,col="red",lty="dotted")
lines(xvec+2,+phiF,lwd=3,col="red",lty="dotted")

dev.copy2pdf(file="phi-periodic-extension.pdf")

# 
