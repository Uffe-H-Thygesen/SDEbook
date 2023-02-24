## Simulation of Geometric Brownian Motion

set.seed(12345)

sigma <- 1

T <- 10
dt <- 0.001

require(SDEtools)

tv <- seq(0,T,dt)

B <- rBM(tv)

pdf("gbm.pdf",width=7,height=4)
par(mfrow=c(1,2))

for(r in c(-0.5,1))
    {
        Y <- B*sigma + (r-0.5*sigma^2)*tv
        X <- exp(Y)
        plot(tv,X,type="l",col="black",xlab=expression(t),ylab=expression(X[t]))
    }

dev.off()
