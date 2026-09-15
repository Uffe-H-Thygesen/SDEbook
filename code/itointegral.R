### Demonstrate problem with int(B*dB), right vs left evaluation
### 

set.seed(1234)
require(SDEtools)

## Fine discretization: Number of steps
N <- 512*16

## Time range
T <- 1
h <- T/N
tv <- seq(0,T,h)

B <- rBM(tv)
dB <- diff(B)

## The Ito integral and the corresponding "right hand evaluation"
Il <- c(0,cumsum(head(B,-1)*dB))
Ir <- c(0,cumsum(tail(B,-1)*dB))

pdf(file="itointegral.pdf",width=7,height=5)
par(mfrow=c(2,1),mar=c(4,4,0.5,2))

plot(tv,B,xlab="",ylab=expression(B[t]),type="l")

plot(c(0,1),range(c(Ir,Il)),type="n",xlab="t",ylab=expression(I[t]))
lines(tv,0.5*B^2,lty="solid",lwd=2)
lines(tv,Ir,lty="solid",col="grey",lwd=2)
lines(tv,Il,lty="solid",col="grey",lwd=2)

### Subsample and do it again
ts <- seq(1,length(B),16)

B <- B[ts]
tv <- tv[ts]
dB <- diff(B)


Il <- c(0,cumsum(head(B,-1)*dB))
Ir <- c(0,cumsum(tail(B,-1)*dB))

lines(tv,Ir,lty="solid",col="black",lwd=1)
lines(tv,Il,lty="solid",col="black",lwd=1)

text(tv[length(tv)],Ir[length(tv)],pos=2,expression(I[t]^R))
text(tv[length(tv)],Il[length(tv)],pos=2,expression(I[t]^L))
text(tv[length(tv)],0.5*B[length(tv)]^2,pos=3,expression(B[t]^2/2))

dev.off()
