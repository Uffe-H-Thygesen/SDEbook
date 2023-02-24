require(SDEtools)

## Function to simulate a Brownian bridge 
BrownianBridge <- function(tv,b=0)
{
    T <- tv[length(tv)]
    B <- rBM(tv)
    BB <- B + (b-B[length(B)])*tv/T
    return(BB)
}

N <- 1000
b <- 1
T <- 2
dt <- 0.01

tv <- seq(0,T,dt)

BB <- replicate(N,BrownianBridge(tv,b))

EBB <- apply(BB,1,mean)
VBB <- apply(BB,1,var)

par(mfrow=c(1,2))

plot(tv,EBB,ylim=c(-1,2))
lines(tv,tv*b/T)
lines(tv,EBB+sqrt(VBB))
lines(tv,EBB-sqrt(VBB))
matplot(tv,BB[,1:2],type="l",lty=1,col=1,add=TRUE)

plot(tv,VBB)
lines(tv,tv*(1-tv/T))


     
