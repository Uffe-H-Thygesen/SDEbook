### Sample path of Ito integral of B w.r.t. itself
### Sample paths as well as 1-sigma limits

T <- 2
h <- 0.001

tvec <- seq(0,T,h)
dB <- rnorm(length(tvec)-1,sd=sqrt(h))
B <- c(0,cumsum(dB))

It <- G*c(0,cumsum(B[-length(tvec)]*dB))

VIt <- 0.5*tvec^2

plot(tvec,It,type="l",ylim=c(-1,1)*2*sqrt(VIt[length(tvec)]),
     xlab="t",ylab=expression(I[t]))
lines(tvec,sqrt(VIt),lty="dashed")
lines(tvec,-sqrt(VIt),lty="dashed")

