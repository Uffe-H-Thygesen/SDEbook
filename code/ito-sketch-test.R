## Sketch illustrating Ito's lemma

graphics.off()
rm(list=ls())

#pdf(file="ito-sketch-test.pdf")

## Define the function h

yoff <- 0.

a <- 4
h <- function(x) x^2 # exp(a*x)
hinv <- function(y) sqrt(y) # 1/a*log(y) # sqrt(y-yoff)
dhdx <- function(x) 2*x # a*exp(a*x) # 2*x

# h <- function(x) yoff+x^2
# hinv <- function(y) sqrt(y-yoff)
# dhdx <- function(x) 2*x

## Initial value of x
x0 <- 0.6

## Range of x coordinates
xmin <- 0
xmax <- 1
xvec <- seq(xmin,xmax,length=1001)

## Drift and diffusion term
f <- function(x) 0
g <- function(x) 1

## Time step
dt <- 0.01

## Density of X after the time step
xnew <- x0+f(x0)*dt
    
fX <- function(x) dnorm(x,mean=xnew,sd=g(x0)*sqrt(dt))

fmax <- fX(xnew)

plot(c(xmin,xmax),c(0,0.5*h(xmax)),type="n",xlab="x",ylab="y")

shade <- 0.5
shading <- rgb(shade,shade,shade,alpha=0.5)

polygon(xvec,2*fX(xvec),col=shading)

plot(h,from=xmin,to=xmax,ylim=c(0,h(xmax)),add=TRUE)


lines(c(x0,x0,xmin),c(0,h(x0),h(x0)),lwd=3)

lines(c(xnew,xnew,xmin),c(0,h(xnew),h(xnew)),lty="dashed")

## Density of Y after the time step 
fY <- function(y) fX(hinv(y))/abs(dhdx(hinv(y)))



yvec <- seq(yoff,h(xmax),length=1001)

fYmax <- max(fY(yvec))

polygon(xmin + 0.25*fY(yvec)/fYmax,yvec,density=10,col=NA)

## fYY <- function(y) dlnorm(y,a*xnew,a*g(x0)*sqrt(dt))
## polygon(xmin + 0.25*fYY(yvec)/fYmax,yvec,density=10,col="red")

X <- rnorm(n=1e5,mean=xnew,sd=sqrt(dt)*g(x0))
Y <- h(X)
EY <- mean(Y)
abline(h=EY,lty="dashed")

# dev.off()
