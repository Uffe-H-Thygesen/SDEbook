## Stability of the double well system dX = (lambda*X-X^3) dt + sigma dB

require(SDEtools)

## Model parameters
sigma <- 0.5
D = 0.5* sigma^2
lambdas <- seq(-1,+1)

## Range for plotting
xmax <- sqrt(max(abs(lambdas)))*2
xmin <- -xmax

## Simulation control
x0 <- 0

T <- 200
dt <- 0.1
tv <- seq(0,T,dt)

B <- rBM(tv)

g <- function(x) sigma

## Simulate the system for each alue of lambda
res <- lapply(lambdas,function(l)
{
    lambda <- l
    f <- function(x) x*(lambda - x^2)
    return(list(lambda=l,sim=euler(f,g,tv,x0,B)))
}
)

pdf("doublewell.pdf",width=6,height=6)

par(mfrow=c(length(lambdas),3),mar=c(4,4,2,2))

## Plot the results for each value of lambda
lapply(res,function(r)
{
    ## The potential 
    U <- function(x) - 0.5 * r$lambda*x^2 + 0.25 * x^4

    ## The drift function
    f <- function(x) r$lambda*x - x^3
    plot(f, from=-2,to=2,ylim=c(-4,4)) # c(-0.25,4))

    ## Add equilibria, coded for their stability
    if(r$lambda<0) points(0,0,pch=16) else points(c(-sqrt(r$lambda),0,sqrt(r$lambda)),c(0,0,0),pch=c(16,1,16))
    grid()
    
    ## Plot sample path
    plot(r$sim$times,r$sim$X,type="l",xlab="t",ylab=expression(X[t]),ylim=c(xmin,xmax))

    ## Plot stationary density, normalized
    phi <- function(x) exp( - U(x)/D )
    Z <- integrate(phi,lower=-Inf,upper=+Inf)$value
    phiN <- function(x) phi(x)/Z
    plot(phiN,from=xmin,to=xmax,xlab="x",ylab=expression(phi))
    
}
)

dev.off()
