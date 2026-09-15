### Figures examining the strong and weak order of the Euler, Mil'shtein and Heun schemes

require(Matrix)
require(SDEtools)

set.seed(123456)

## Number of replicates
Nsim <- 1e5

## Parameters for geometric Brownian motion
r <- 0.5
sigma <- 0.5
x <- 1

## Simulation interval and finest partition
T <- 1
NstepsMax <- 2^10
hsteps <- 6
tv <- seq(0,T,length=NstepsMax+1)

hmin <- T/NstepsMax

## Increments in Brownian motion - each column is a realization, each row a time step
dB <- array(sqrt(hmin)*rnorm(Nsim*NstepsMax),c(NstepsMax,Nsim))

## Analytical solution of GBM at time T
Xa <- x*exp((r-0.5*sigma^2)*T+sigma*apply(dB,2,sum))

## Array for the approximated solutions. One row for each grid size;
## one column for each realization
Xm <- Xh <- Xe <- array(0,c(hsteps,Nsim))

## Ito and Stratonovich Drift and noise for the Milshtein and  Heun method
f <- function(x) x*r
fS <- function(x) x*(r-0.5*sigma^2)

g <- function(x) Diagonal(Nsim,sigma*as.numeric(x))
gp <- function(x) Diagonal(Nsim,sigma)

milshtein <- function()
{
    X <- rep(x,Nsim)
    for(i in 2:length(tv))
        X <- X + f(X)*h + g(X) %*% as.numeric(dB[i-1,]) + 0.5*gp(X) %*% g(X) %*% ( as.numeric(dB[i-1,]^2 - h))
    return(X)
}

h <- hmin
hv <- hmin*2^(0:(hsteps-1))

for(i in 1:hsteps)
  {
    print(i)

    B <- apply(dB,2,function(db)c(0,cumsum(db)))

    ## Determine discretized solution using the
    ## "analytical" solution of the discretized equation
    Xe[i,] <- x*apply(1+r*h+sigma*dB,2,prod)

    ## Mil'shtein method
    Xm[i,] <- x*apply(1+r*h+sigma*dB+0.5*sigma^2*(dB^2-h),2,prod)

    ## Heun method
    ## Xh[i,] <- heun(f=fS,g=g,times=tv,x0=rep(x,Nsim),B=B)$X[length(tv),]
    dum <- (r-0.5*sigma^2)*h+sigma*dB
    Xh[i,] <- x*apply(1+dum+0.5*dum^2,2,prod)
    
    ## Sub-sample time and Brownian motion to prepare for coarser simulation 
    h <- 2*h
    dB <- apply(B[seq(1,dim(B)[1],2),],2,diff)
    tv <- tv[seq(1,length(tv),2)]
  }

## Compute error, for each time step and each realization
Xeerr <- Xe - outer(rep(1,hsteps),Xa)
Xherr <- Xh - outer(rep(1,hsteps),Xa)
Xmerr <- Xm - outer(rep(1,hsteps),Xa)

## Compute strong error - mean abs error - for each time step
meErr <- apply(abs(Xeerr),1,mean)
mhErr <- apply(abs(Xherr),1,mean)
mmErr <- apply(abs(Xmerr),1,mean)

### Weak error.

## Test function
hfun <- function(x) x

## Compute weak error. Dividing into batches for smaller statistical error
nBatch <- 10
Batch <- rep(1:nBatch,rep(Nsim/nBatch,nBatch))
weakerror <- function(x) tapply(hfun(x)-hfun(Xa),Batch,mean)

## Compute weak error for each batch and each time step 
wXeerr <- apply(Xe,1,weakerror) 
wXmerr <- apply(Xm,1,weakerror) 
wXherr <- apply(Xh,1,weakerror) 

## Mean weak error for each time step
mwXeerr <- apply(wXeerr,2,function(x)mean(abs(x)))
mwXmerr <- apply(wXmerr,2,function(x)mean(abs(x)))
mwXherr <- apply(wXherr,2,function(x)mean(abs(x)))

## Richardson extrapolation; order 2, single step
R2 <- function(e) (2*head(e,-1)-tail(e,-1))

## Perform  Richardson extrapolation in each batch 
R2Xeerr <- t(abs(apply(wXeerr,1,R2)))
R2Xmerr <- t(abs(apply(wXmerr,1,R2)))
R2Xherr <- t(abs(apply(wXherr,1,R2)))

## Mean Richardson error for each time step, averaging over batches
mR2Xeerr <- apply(R2Xeerr,2,function(x)mean(abs(x)))
mR2Xmerr <- apply(R2Xmerr,2,function(x)mean(abs(x)))
mR2Xherr <- apply(R2Xherr,2,function(x)mean(abs(x)))

## PLOTS 

pdf(file="figStrongError.pdf", width=5,height=4)

plot(hv,meErr,log="xy",xlab="Time step h",ylab="Strong error",
     ylim=range(c(meErr,mhErr,mmErr)),pch=3)
points(hmin*2^(0:(hsteps-1)),mmErr,pch=4)
points(hmin*2^(0:(hsteps-1)),mhErr,pch=1)

# Add line with slope 0.5
h1 <- 2*hmin
h2 <- hmin*2^(hsteps-2)

e1 <- meErr[1]*2
e2 <- sqrt(h2/h1)*e1

e3 <- (mmErr[2]+mhErr[2])/2
e4 <- e3*h2/h1

lines(c(h1,h2),c(e1,e2),lty=1)
lines(c(h1,h2),c(e3,e4),lty=2)
legend("bottomright",
       c("E-M","Milshtein","Heun"),
       lty=c(-1,-1,-1),pch=c(3,4,1))

dev.off()


ylim <- range(c(mR2Xeerr,mwXeerr,mR2Xmerr,mwXmerr,mR2Xherr,mwXherr))

pdf(file="figWeakError.pdf",width=7,height=5)

plot(hv,mwXeerr,
     log="xy",xlab="Time step h",ylab="Weak error",ylim=ylim,pch=3)
points(hv,mwXmerr,pch=4)
points(hv,mwXherr,pch=1)

## Add Richardson extrapolated errors
points(head(hv,-1),mR2Xherr,pch=16)


e1 <- abs(mwXeerr[1])*4
e2 <- h2/h1*e1

lines(c(h1,h2),c(e1,e2))

legend("bottomright",
       c("E-M","Milshtein","Heun","Richardson"),
       lty=c(-1,-1,-1,-1),pch=c(3,4,1,16))

dev.off()
