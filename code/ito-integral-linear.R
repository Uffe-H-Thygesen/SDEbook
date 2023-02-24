### Sample path of Ito integral associated with linear system
### Sample paths as well as 1-sigma limits

pdf("ito-integral-linear.pdf")

A <- -1
G <- 1

T <- 2
h <- 0.001

tvec <- seq(0,T,h)
dB <- rnorm(length(tvec)-1,sd=sqrt(h))
B <- c(0,cumsum(dB))

iexpmAsGdB <- G*c(0,cumsum(exp(-A*tvec[-length(tvec)])*dB))

Xt <- exp(A*tvec)*iexpmAsGdB


iexpm2AsG2ds <- G^2*c(0,cumsum(exp(-2*A*tvec[-length(tvec)])*h))

VXt <- exp(2*A*tvec) * iexpm2AsG2ds

plot(tvec,Xt,type="l",ylim=c(-1,1)*2*sqrt(VXt[length(tvec)]),
     xlab="t",ylab=expression(X[t]))
lines(tvec,sqrt(VXt),lty="dashed")
lines(tvec,-sqrt(VXt),lty="dashed")

dev.off()
