## Sketch illustrating Ito's lemma

graphics.off()
rm(list=ls())

pdf(file="ito-sketch.pdf")

## Define the function h
a <- 4
h <- function(x) exp(a*x)
hinv <- function(y) 1/a*log(y) # sqrt(y-yoff)
dhdx <- function(x) a*exp(a*x) # 2*x

## Initial value of x
x0 <- 0.6

## Range of x coordinates
xmin <- 0
xmax <- 1
xvec <- seq(xmin,xmax,length=1001)

## Density of X 
mu <- 0.6
sigma <- 0.1
fX <- function(x) dnorm(x,mean=mu,sd=sigma)

## Density of Y 
fY <- function(y) fX(hinv(y))/abs(dhdx(hinv(y)))

## Compute E(Y) in two ways
EY.1 <- integrate(function(y)fY(y)*y,lower=0,upper=Inf)$value
EY.2 <- h(mu)*exp(0.5*a^2*sigma^2)

## Setup the plot
plot(c(xmin,xmax),c(0,0.5*h(xmax)),type="n",xlab="x",ylab="y")

shade <- 0.5
shading <- rgb(shade,shade,shade,alpha=0.5)

## Add density of X (scaled)
polygon(xvec,2*fX(xvec),col=shading)

## Add graph of function h
plot(h,from=xmin,to=xmax,ylim=c(0,h(xmax)),add=TRUE)

## Add lines to the medians
lines(c(mu,mu,xmin),c(0,h(mu),h(mu)),lwd=3)

## Compute the maximum pdf of Y in two ways
yvec <- seq(h(xmin),h(xmax),length=1001)
fYmax.1 <- max(fY(yvec))
Ymode <- exp(4*mu-0.5*a^2*sigma^2)
fYmax.2 <- fY(Ymode)

## Add density of Y (scaled)
polygon(xmin + 0.25*fY(yvec)/fYmax.2,yvec,density=10,col=NA)

## Add line for mean of Y
abline(h=EY.1,lty="dashed")

dev.off()
