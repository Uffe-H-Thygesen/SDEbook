### Illustration for "transport" chapter.

set.seed(123456)


### Figure for Gauss curves

pdf(file="figGauss.pdf",width=7,height=4)

D <- 1
tvec <- c(1,5,10)

xmin <- -10
xmax <- -xmin

ltys <- c("solid","dashed","dotted")

plot(function(x)dnorm(x,sd=sqrt(2*D*tvec[1])),from=xmin,to=xmax,
     xlab="Position [m]",lty=ltys[1],ylab="Concentration [1/m]")

for(i in 2:length(tvec))
  plot(function(x)dnorm(x,sd=sqrt(2*D*tvec[i])),lty=ltys[i],
       from=xmin,to=xmax,add=TRUE)


legend(4.5,0.25,legend=paste("Time =",tvec,"s"),lty=ltys)

dev.off()



### Simulation of random walk and Brownian motion

k <- 1
h <- 1
p <- 0.25

N <- 10000

U <- runif(N)

dX <- numeric(N) - (U<p) + (U>(1-p))

X <- c(0,cumsum(dX))
T <- (0:N)*h

Nplot1 <- 100


pdf("figRW.pdf",height=4)
par(mfrow=c(1,2))

plot(T[1:Nplot1],X[1:Nplot1],type="s",xlab="Time [ns]",ylab=expression(paste("Position [",mu,"m]")))
plot(T,X,type="s",xlab="Time [ns]",ylab=expression(paste("Position [",mu,"m]")))
dev.off()

# dev.copy2pdf(file="figRW.pdf")
