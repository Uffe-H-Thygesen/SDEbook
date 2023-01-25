## Simulation of Geometric Brownian Motion

set.seed(12345)

r <- 1
sigma <- 2
N <- 5

T <- 10
dt <- 0.001

require(SDEtools)
tv <- seq(0,T,dt)

B <- rvBM(tv,N)
Y <- B*sigma + (r-0.5*sigma^2)*tv
X <- exp(Y)

pdf("gbm-sim.pdf",width=6,height=3.5)
matplot(tv,X,type="l",col="black",xlab=expression(t),ylab=expression(X[t]),ylim=c(0,20))
dev.off()
