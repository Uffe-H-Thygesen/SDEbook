graphics.off()
rm(list=ls())

require(SDEtools)

k <- 1
c <- 0.1
sigma <- 1
s <- 1

h <- 0.05
T <- 100

tv <- seq(0,T,h)

B <- rBM(tv)
W <- rBM(tv)

A <- array(c(0,-k,1,-c),c(2,2))
G <- array(c(0,sigma),c(2,1))
C <- array(c(1,0),c(1,2))
D <- s/sqrt(h)
DD <- D^2

sys <- dLinSDE(A,G,h)
## sys$eAt <- as.matrix(sys$eAt)

F <- chol(sys$St)

Xsim <- array(0,c(2,length(tv)))

for(i in 2:length(tv))
{
    Xsim[,i] <- sys$eAt %*% Xsim[,i-1] + t(F) %*% rnorm(2)
}

Ysim <- c(0,Xsim[1,-1] + s/h*diff(W))

par(mfrow=c(2,1))
plot(tv,Xsim[1,],type="l")
points(tv,Ysim)

plot(tv,Xsim[2,],type="l")

## Filtering
S <- array(0,c(2,2,length(tv)))
Xhat <- array(0,c(2,length(tv)))
           
for(i in 2:length(tv))
{
    ## Time update
    Xhat[,i] <- sys$eAt %*% Xhat[,i-1]
    S[,,i] <- sys$eAt %*% S[,,i-1] %*% t(sys$eAt) + sys$St

    ## Kalman gain
    K <- S[,,i] %*% t(solve(C %*% S[,,i] %*% t(C) + DD,C) )

    ## Data update
    Xhat[,i] <- Xhat[,i] + K %*% (Ysim[i] - C %*% Xhat[,i])

    S[,,i] <- S[,,i]  - K %*% C %*% S[,,i]
}

dev.new()
my.conf.int <- function(t,xmean,sd)
{
    polygon(c(t,rev(t)),c(xmean+sd,rev(xmean-sd)),col="gray")
}

par(mfrow=c(2,1))

sdX <- sqrt(S[1,1,])
xlim <- range(c(Xsim[1,],Xhat[1,]+sdX,Xhat[1,]-sdX))
plot(range(tv),range(xlim),type="n")
my.conf.int(tv,Xhat[1,],sdX)
lines(tv,Xhat[1,])
points(tv,Xsim[1,])


sdV <- sqrt(S[2,2,])
xlim <- range(c(Xsim[2,],Xhat[2,]+sdV,Xhat[2,]-sdV))
plot(range(tv),range(xlim),type="n")
my.conf.int(tv,Xhat[2,],sdV)
lines(tv,Xhat[2,])
points(tv,Xsim[2,])

print(var(t(Xsim)))
print(var(t(Xsim-Xhat)))
print(S[,,length(tv)])

