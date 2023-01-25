### Figure for advective-diffusive transport

C0 <- function(x) dnorm(x,mean=0,sd=1)
CT <- function(x) dnorm(x,mean=10,sd=3)


pdf("advdiff-fig.pdf",width=7,height=4)

myplot <- function(f,from,to,col,...)
  {
    plot(f,from=from,to=to,...)
    xx <- seq(from,to,length=1000)
    ff <- f(xx)
    polygon(c(xx,to,from),c(ff,0,0),col=col)
  }

myplot(C0,from=-3,to=20,col="darkgray",ylim=c(0,0.4),xlab="x",ylab="C")
myplot(CT,from=-3,to=20,col="lightgray",ylim=c(0,0.4),add=TRUE)

text(2.5,0.2,"C(x,0)")
text(17,0.05,"C(x,T)")

dev.off()
