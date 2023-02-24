### Demonstrate problem with int(B*dB), right vs left evaluation
### 
###

set.seed(1234)

## Fine discretization: Number of steps
N <- 512*16

## Time range
T <- 1

h <- T/N

dB <- rnorm(N,sd=sqrt(h))
B <- c(0,cumsum(dB))

Il <- c(0,cumsum(B[-(N+1)]*dB))
Ir <- c(0,cumsum(B[-1]*dB))

tvec <- seq(0,T,length=N+1)

pdf(file="itointegral.pdf",width=7,height=5)


par(mfrow=c(2,1),mar=c(4,4,0.5,2))

plot(tvec,B,xlab="",ylab=expression(B[t]),type="l")

plot(c(0,1),range(c(Ir,Il)),type="n",xlab="t",ylab=expression(I[t]))

lines(tvec,0.5*B^2,lty="solid",lwd=2)
lines(tvec,Ir,lty="solid",col="grey",lwd=2)
lines(tvec,Il,lty="solid",col="grey",lwd=2)

### Subsample
ts <- seq(1,length(B),16)

B <- B[ts]
tvec <- tvec[ts]
dB <- diff(B)


Il <- c(0,cumsum(B[-length(B)]*dB))
Ir <- c(0,cumsum(B[-1]*dB))

lines(tvec,Ir,lty="solid",col="black",lwd=1)
lines(tvec,Il,lty="solid",col="black",lwd=1)
# lines(tvec,0.5*B^2,lty="solid",lwd=1)

text(tvec[length(tvec)],Ir[length(tvec)],pos=2,expression(I[t]^R))
text(tvec[length(tvec)],Il[length(tvec)],pos=2,expression(I[t]^L))
text(tvec[length(tvec)],0.5*B[length(tvec)]^2,pos=3,expression(B[t]^2/2))

dev.off()
