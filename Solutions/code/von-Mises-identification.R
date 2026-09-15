## Simulation / re-estimation of parameters using Girsanov

require(SDEtools)

## Simulation model
a <- 1
theta <- 0
sigma <- 1

f <- function(x) a*sin(theta-x)
g <- function(x) sigma
    
T <- 100
dt <- 0.1
times <- seq(0,T,dt)
    
x0 <- 0

## Generate data set
X <- euler(f,g,times,x0)$X

## Function to compute the profile likelihood of theta as well as
## the ML estimate of a, given theta
##
## With Girsanov's theorem, the likelihood of a pair (theta,a) is
## 
profile.likelihood <- function(th)
{
    idX <- tail(itointegral(sin(th-X)/g(X)^2,X),1)
    idt <- tail(itointegral(sin(th-X)^2/g(X)^2,times),1)
    return(c(logL = 0.5*idX^2/idt,ahat=idX/idt))
}

## Tabulate the profile likelihood
thetas <- seq(-pi,pi,length=1001)
profL <- sapply(thetas,profile.likelihood)

## Identify the Maximum Likelihood estimate of theta
## (Note: We could refine this estimate by using it as a starting guess for
## optim)
ihat <- which.max(profL[1,])
thetahat <- thetas[ihat]
ahat <- profL[2,ihat]

## Use symmetry to get estimates with a>0 and -pi <= theta < pi
if(ahat <0) {
    ahat <- -ahat
    thetahat <- thetahat + pi
}
thetahat <- (thetahat + pi) %% (2*pi) - pi

## Assess the variance of the theta-estimate from the curvature of log L
n <- length(thetas)
dth <- diff(thetas[1:2])

curvature <- (profL[1,(ihat-2)%%n + 1] -2*profL[1,ihat]  + profL[1,(ihat %% n)+1] ) / dth^2
Vthetahat <- -1/curvature

## Plot data, profile likelihood, and conditional estimate of parameter "a"
## (Note the symmetry)
##
## Add estimates and confidence interval
par(mfrow=c(3,1))
plot(times,X,type="l",xlab="Time",ylab="X",main="Data")
plot(thetas,profL[1,],type="l",xlab=expression(theta),
     ylab="log(L)",main="Profile likelihood of theta")

points(thetahat,max(profL[1,]),pch=16,cex=3)
abline(v=thetahat-2*sqrt(Vthetahat),lty="dashed")       
abline(v=thetahat+2*sqrt(Vthetahat),lty="dashed")
abline(h=max(profL[1,])-2,lty="dashed")

plot(thetas,profL[2,],type="l",xlab=expression(theta),
     ylab=expression(hat(a)),main="Most likely 'a' coefficient")
points(thetahat,ahat,cex=3,pch=16)

