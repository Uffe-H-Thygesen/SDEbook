## Simulation of Brownian Motion with drift

set.seed(12345)

sigma <- 1
u <- 1

T <- 100
dt <- 0.001
T1 <- 0.01

require(SDEtools)

tv <- c(seq(0,T1,length=1001),seq(T1,T,length=1001))

B <- rBM(tv)

pdf("bm-drift.pdf",width=7,height=4)
par(mfrow=c(1,2))

X <- B*sigma + tv*u

ylim <- range(X)

I <- (tv<=T1)

plot(tv[I],X[I],type="l",col="black",xlab=expression(t),ylab=expression(X[t]))
plot(tv,X,type="l",col="black",xlab=expression(t),ylab=expression(X[t]))

dev.off()
