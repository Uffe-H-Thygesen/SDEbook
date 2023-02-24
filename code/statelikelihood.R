graphics.off()
rm(list=ls())

pdf(file="statelikelihood.pdf",width=5,height=2.5)

par(mfrow=c(1,3))

xmin <- 0
xmax <- 10

ylab <- expression(l[i](x))

l <- function(x) dnorm(3,mean=x,sd=sqrt(0.5))

plot(l,from=xmin,to=xmax,ylab=ylab,main="Gaussian")

print(integrate(l,lower=0,upper=10))

l <- function(x) dpois(3,lambda=x)
plot(l,from=xmin,to=xmax,ylab=ylab,main="Poisson")

print(integrate(l,lower=0,upper=10))

l <- function(x) as.numeric(abs(x-3)<0.5)
plot(l, ,from=xmin,to=xmax,ylab=ylab,main="Round-off")

print(integrate(l,lower=0,upper=10))

dev.off()
