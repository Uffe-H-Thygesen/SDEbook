## Modified Euler for Geometric Brownian Motion

## Parameters
r <- 1
sigma <- 0.25
T <- 1
x0 <- 1


stochrk <- function (f, g, times, x0, B = NULL, p = function(x) x) 
{
    nx <- length(x0)
    nt <- length(times)
    if (length(formals(f)) == 1) {
        ff <- function(t, x) f(x)
    }
    else {
        ff <- f
    }
    if (length(formals(g)) == 1) {
        ggg <- function(t, x) g(x)
    }
    else {
        ggg <- g
    }
    g0 <- ggg(times[1], x0)
    if (is.null(dim(g0))) {
        nB <- 1
        gg <- function(t, x) matrix(ggg(t, x), nrow = nx, ncol = 1)
    }
    else {
        nB <- ncol(g0)
        gg <- ggg
    }
    if (is.null(B)) {
        B <- rvBM(times, nB)
    }
    if (!is.matrix(B)) 
        B <- matrix(B, nrow = nt, ncol = nB)
    dB <- apply(B, 2, diff)
    X <- array(NA, c(nt, nx))
    X[1, ] <- x0
    dt <- diff(times)

    ## RK 4
    s <- 4
    cv <- c(0,0.5,0.5,1)
    bv <- c(1,2,2,1)/6
    am <- matrix(c(0 , 0 , 0 ,0 ,
                   0.5,0 , 0 ,0 ,
                   0,  0.5,0 ,0 ,
                   0,  0,  1, 0 ), byrow=TRUE,nrow=4)

    ## Heun
    ## s <- 2
    ## cv <- c(0,1)
    ## bv <- c(0.5,0.5)
    ## am <- matrix(c(0,0,
    ##                1,0), byrow=TRUE,nrow=2)
    
    k <- array(0,c(s,nx))

    for (i in 1:(nt - 1)) {
        fg <- function(t,x) as.numeric(ff(t,x) + gg(t,X[i,]) %*% as.numeric(dB[i,]) / dt[i])
        ## fg <- function(t,x) as.numeric(ff(t,x) + gg(t,x) %*% as.numeric(dB[i,]) / dt[i])

        for(j in 1:s)
        {
            pred <- p(X[i,] +  as.numeric( am[j,] %*% k * dt[i] ))
            k[j,] <- as.numeric(fg(times[i]+cv[j]*dt[i],pred))
        }
        
        X[i + 1, ] <- as.numeric(p(X[i, ] + bv %*% k * dt[i]))
    }
    colnames(X) <- names(x0)
    return(list(times = times, X = X))
}


## Number of realizations
Nsim <- 1000

## Finest step size
hmin <- 2^(-12)
nh <- 8

hvec <- hmin * 2^seq(0,nh-1,1)

require(SDEtools)

tv <- seq(0,T,length=T/hmin + 1)
B <- rvBM(tv,Nsim)

dB <- apply(B,2,diff)

Xa <- x*exp((r-0.5*sigma^2)*T+sigma*B[length(tv),])

rk <- rkerr <- EM <- Mod <- ModErr <- EMerr <- array(NA,c(length(hvec),Nsim))

h <- hmin 

f <- function(x) r*x
fS <- function(x) (r-0.5*sigma^2)*x
g <- function(x) sigma*x

for(i in 1:nh)
{
    Xe <- x*apply(1+r*h+sigma*dB,2,prod)
    Xm <- x*apply(exp(r*h)+sigma*dB,2,prod)

    EM[i,] <- Xe
    Mod[i,] <- Xm
    rk[i,] <- stochrk(f, g, tv, rep(x,1))$X[length(tv),]
    
    EMerr[i,] <- Xe - Xa
    ModErr[i,] <- Xm - Xa
    rkerr[i,] <- rk[i,] - Xa
    
    B <- B[seq(1,nrow(B),2),]
    dB <- apply(B,2,diff)
    h <- 2*h
}

## Compute strong error
meanabsEMerr <- apply(EMerr,1,function(x)mean(abs(x)))
meanabsModErr <- apply(ModErr,1,function(x)mean(abs(x)))
meanabsrkErr <- apply(rkerr,1,function(x)mean(abs(x)))

## Compute weak error
## Test function
k <- function(x) x^2
weakEMerr <- apply(k(EM),1,mean) - mean(k(Xa))
weakModErr <- apply(k(Mod),1,mean) - mean(k(Xa))
weakrkErr <- apply(k(rk),1,mean) - mean(k(Xa))

## Richardson extrapolation
R2EMerr  <- (2*head(weakEMerr ,-1)-tail(weakEMerr ,-1))/(2-1)
R2ModErr <- (2*head(weakModErr,-1)-tail(weakModErr,-1))/(2-1)
R2ModErr <- (2*head(weakModErr,-1)-tail(weakModErr,-1))/(2-1)

meanabsR2EMerr <- apply(R2EMerr,1,function(x)mean(abs(x)))
meanabsR2ModErr <- apply(R2ModErr,1,function(x)mean(abs(x)))

## Plots 
ylim <- range(c(meanabsEMerr,meanabsModErr,meanabsR2EMerr,meanabsR2ModErr))

plot(hvec,meanabsEMerr,log="xy",ylim=ylim)
points(hvec,meanabsModErr,pch=2)
points(head(hvec,-1),meanabsR2EMerr,pch=3,col="red")
points(head(hvec,-1),meanabsR2ModErr,pch=4,col="red")

hs <- c(hvec[2],hvec[nh-1])
es <- min(meanabsModErr)*sqrt(hs/hmin)/sqrt(2)
lines(hs,es,lwd=2)

