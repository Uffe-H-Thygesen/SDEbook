require(SDEtools)
require(fields)
require(MASS)

## Transition probabilities for the double well model

## Potential and noise
gamma <- 10
sigma <- 1

u <- function(x) gamma*x*(1-x^2)
D <- function(x) 0.5*sigma^2

## Spatial grid
Xm <- 2
ni <- 201

xi <- seq(-Xm,Xm,length=ni)
xc <- cell.centers(xi)

dx <- xi[2] - xi[1]

## Generator
G <- fvade(u,D,xi,'r')

## Stationary distribution (nornmalized to sum to one)
phi <- StationaryDistribution(G)

## Time grid
tv <- 10^seq(-4,3,length=51)

## Inital density: Atom at x0
x0 <- 0.25
phi0 <- 0*xc
phi0[ sum(xc<x0) ] <- 1

## Convert probabilities to credibilities
pdf2cred <- function(pdf)
{
    I <- order(pdf)
    cred <- pdf
    cred[I]  <- cumsum(pdf[I])
    return(cred)
}

## Compute forward probabilities
PHI <- as.matrix(do.call(rbind,sapply(tv,function(t)phi0 %*% expm(G*t))))
CRED <- t(apply(PHI,1,pdf2cred))
CDF <- t(apply(PHI,1,cumsum))

## Transition probabilility matrix
P <- as.matrix(expm(G*0.1))

## Backward likelihood of ending near x=0.75
kT <- 0*phi0
kT[ sum(xc<0.75) ] <- 1/dx

PSI <- as.matrix(do.call(cbind,sapply(tv,function(t) expm(G*t) %*% kT)))

PSIr <- PSI[,dim(PSI)[2]:1]

logPSI <- log(abs(PSIr))
logPSI[logPSI < -3 ] <- -3


## Solution of ODE (drift only) for comparison
model <- function(t,x,p) return(list(gamma*x*(1-x^2)))
require(deSolve)
x0 <- xi
xt <- sapply(xi,function(x)(ode(x,c(0,0.1),model,NULL))[2,2])

## Graphics

col <- NULL #gray.colors(32)
width = 5

my.image <- function(x,y,z,col=col,xlab,ylab,main,xl=NULL,yl=NULL)
{
    layout(matrix(1:2,nrow=1),width=c(0.75,0.25))
    par(mar=c(5,4,4,0)+0.1)
    plot(range(x),range(y),type="n",xlab=xlab,ylab=ylab,main=main,xaxs="i",yaxs="i")
    rangex <- range(x)
    rangey <- range(y)
    rangez <- range(z)
    scale <- function(z) (z-rangez[1])/diff(rangez)
    rasterImage(t(scale(z[,ncol(z):1])),rangex[1],rangey[1],rangex[2],rangey[2])

    if(!is.null(xl)) lines(xl,yl)
    plot(c(0,1),rangez,type="n",xaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
    rasterImage(seq(1,0,length=101),0,rangez[1],1,rangez[2])
}


pdf(file="fwdKolmo.pdf",width=width)
my.image(log10(tv),xc,CRED,col=col,xlab=expression(log[10]*t),ylab="y",main=expression(pi(t,y)))
dev.off()

pdf(file="bwdKolmo.pdf",width=width)
my.image(log10(tv),xc,t(logPSI),col=col,
      xlab=expression(-log[10]*(t-s)),ylab="x",main=expression(psi(s,x)))
tvlab <- seq(-4,2,length=7)
axis(1,tvlab,labels=rev(tvlab))
dev.off()

require(latex2exp)

pdf(file="transprob.pdf",width=width)
my.image(xc,xc,(P)/dx, col=col,
         xlab="x",ylab="y",
         main=TeX("$p(0\\rightarrow t,x\\rightarrow y)$"),
         xl=x0,
         yl=xt)


dev.off()

