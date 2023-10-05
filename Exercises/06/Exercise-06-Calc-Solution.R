## -----------------------------------------------------------------------------
## Basic SDE
f <- function(x) x*(1-x)
g <- function(x) sigma*x

## Transform, its derivatives, and its inverse
h <- function(x) log(x)/sigma
dhdx <- function(x) 1/x/sigma
dh2dx2 <- function(x) -1/x^2/sigma

hi <- function(y) exp(y*sigma)


## -----------------------------------------------------------------------------
## Use package SDEtools for rBM, rvBM, stochint, and covariation
require(SDEtools) 

## Helper for simulating X Y and Z, given an initial condition, a time grid, and 
## possibly a sample path of the Brownian motion 
simXYZ <- function(X0,tv,B=NULL)
{
    ## If no sample path of B is given, simulate one
    if(is.null(B)) B <- rBM(tv)

    nt <- length(tv)
    dt <- diff(tv)

    dB <- diff(B)
    
    Y <- X <- numeric(nt)
    X[1] <- x0
    
    Y[1] <- h(X[1])
    
    for(i in 1:(nt-1))
    {
        dX <- f(X[i])*dt[i] + g(X[i]) * dB[i]
        X[i+1] <- X[i] + dX

        ## Use the general expression:
        ##
        ## Y[i+1] <- Y[i] + dhdx(X[i])*dX + 0.5*dh2dx2(X[i])*g(X[i])*g(X[i])*dt[i]
        ##
        ## or the specific one: 

        Y[i+1] <- Y[i] + ((1-exp(Y[i]*sigma))/sigma  - sigma/2 ) * dt[i] + dB[i]
    }
    
    ## Compute Z by stochastic integration
    Z <- Y[1] +  stochint(dhdx(X),X) + 0.5*stochint(dh2dx2(X),QuadraticVariation(X))

    return(list(X=X,Y=Y,Z=Z))
}


## -----------------------------------------------------------------------------
## Compare different algorithms
sigma <- 0.2
x0 <- 0.1

tv <- seq(0,10,0.01)

sim <- simXYZ(X0,tv)

plot(tv,sim$X,type="l",
     ylim=c(0,max(c(sim$X,sim$Y,sim$Z))),
     xlab="Time",
     ylab="Abundance")
lines(tv,hi(sim$Y),col="red")
lines(tv,hi(sim$Z),col="blue")


## -----------------------------------------------------------------------------
print(max(abs(sim$X-hi(sim$Y))))
print(max(abs(sim$X-hi(sim$Z))))


## -----------------------------------------------------------------------------
## The effect of sigma on the stationary distribution
sigmas <- seq(0,2,0.25)
vars <- means <- numeric(length(sigmas))

## Long simulation
tv <- seq(0,100,0.01)

## Use the same Brownian motion 
B = rBM(tv)

par(mfrow=c(3,3))

for(j in 1:length(sigmas))
{
    sigma <- sigmas[j]

    sim <- simXYZ(X0,tv,B)

    ## Base statistics only on the tail 
    Xtail <- tail(sim$X,round(0.9*length(tv)))
    means[j] <- mean(Xtail)
    vars[j] <- var(Xtail)

    ## Plot the histogram
    ## Fix the break points so it is easier to compare between different values of sigma
    Xmax <- 2.5 
    hist(pmin(Xtail,Xmax),breaks=seq(0,Xmax,0.1),main=sigma,xlab="x",freq=FALSE)
}


## -----------------------------------------------------------------------------
plot(sigmas,means,ylim=c(0,1.5),xlab=expression(sigma),ylab="Mean +/- std.dev")
lines(sigmas,means+sqrt(vars),lty="dashed")
lines(sigmas,means-sqrt(vars),lty="dashed")


## -----------------------------------------------------------------------------
A <- array(c(0,-1,1,0),c(2,2))
sigma <- 1
G <- sigma*diag(c(1,1))

tvec <- seq(0,10,0.001)

B <- rvBM(tvec,n=2)
sim <- euler(f=function(x) A %*% x,g=function(x) G,times=tvec,x0=numeric(2),B=B)

S <- apply(sim$X^2,1,sum)
S1 <- 2*stochint(sigma*sim$X[,1],B[,1]) + 2*stochint(sigma*sim$X[,2],B[,2]) + 2*sigma^2*tvec

plot(tvec,S,type="l")

lines(tvec,S1,col="red")

