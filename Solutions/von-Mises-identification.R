## Simulation / re-estimation of parameters using Girsanov

require(SDEtools)

asdf <- function(){
    a <- 1
    theta <- 0
    sigma <- 1.

    f <- function(x) a*sin(theta-x)
    g <- function(x) sigma
    
    T <- 100
    dt <- 0.1
    times <- seq(0,T,dt)
    
x0 <- 0

X <- euler(f,g,times,x0)$X

profile.likelihood <- function(th)
{
    idX <- tail(itointegral(sin(th-X)/g(X)^2,X),1)
    idt <- tail(itointegral(sin(th-X)^2/g(X)^2,times),1)
    return(c(logL = 0.5*idX^2/idt,ahat=idX/idt))
}

thetas <- seq(0,pi,length=1001)

profL <- sapply(thetas,profile.likelihood)

## par(mfrow=c(3,1))
## plot(times,X,type="l")
## plot(thetas,profL[1,],type="l")
## plot(thetas,profL[2,],type="l")

ihat <- which.max(profL[1,])
thetahat <- thetas[ihat]
ahat <- profL[2,ihat]

if(ahat <0) {
    ahat <- -ahat
    thetahat <- thetahat + pi
}

thetahat <- (thetahat + pi) %% (2*pi) - pi

n <- length(thetas)
dth <- diff(thetas[1:2])

curvature <- (profL[1,(ihat-2)%%n + 1] -2*profL[1,ihat]  + profL[1,(ihat %% n)+1] ) / dth^2
Vthetahat <- -1/curvature

return(c(thetahat,Vthetahat))
}


qwer <- replicate(100,asdf())

mean(qwer[1,])
var(qwer[1,])
mean(qwer[2,])
