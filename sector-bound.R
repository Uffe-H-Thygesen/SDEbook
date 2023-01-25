C <- 1

f <- function(x) -(x+1)*exp(x)/(1+exp(x))

g <- function(x) 2*tanh(x)


Xmax <- 4

bound <- function(x) C*(1+abs(x))

plot(Xmax*c(-1,1),c(-1,1)*bound(Xmax),type="n",xlab="x",ylab=expression(f(x)*", "*g(x)))

Xvec <- c(-Xmax,0,Xmax,Xmax,-Xmax)
Yvec <- c(bound(c(-Xmax,0,Xmax,Xmax,Xmax)))

polygon(Xvec,Yvec,border=NA,col="grey")
polygon(Xvec,-Yvec,border=NA,col="grey")

plot(f,from=-Xmax,to=Xmax,add=TRUE,lwd=2,lty=1)
plot(g,from=-Xmax,to=Xmax,add=TRUE,lwd=2,lty=3)

text(-Xmax,f(-Xmax),pos=1,expression(f(x)))
text(-Xmax,g(-Xmax),pos=1,expression(g(x)))

grid()

dev.copy2pdf(file="sector-bound.pdf")
