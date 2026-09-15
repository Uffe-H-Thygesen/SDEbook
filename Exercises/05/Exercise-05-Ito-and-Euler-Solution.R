## -----------------------------------------------------------------------------
itointegral <- function(G,B) c(0,cumsum(head(G,-1)*diff(B)))


## -----------------------------------------------------------------------------
## Time grid
tvec <- seq(0,2*pi,length=51)

## Compute integral
I <- itointegral(cos(tvec),sin(tvec))

## Plots
plot(tvec,I,type="l",col="red",
     xlab="Time",ylab="Integral of cos w.r.t. sin")
plot(function(t)t/2+sin(2*t)/4,col="black",from=0,to=2*pi,add=TRUE)

legend(x="topleft",leg=c("Left","Analytical"),
       lty="solid",col=c("red","black"))


## -----------------------------------------------------------------------------
## Time grid
tvec <- seq(0,100,0.5) ## Use a finer time step for increased accuracy

## Simulate Brownian motion
rBM <- function (times) cumsum(rnorm(length(times),sd=sqrt(diff(c(0,times)))))
Bvec <- rBM(tvec)

## Compute integral
I <- itointegral(Bvec,Bvec)

## Compute also "right hand rule"
rightintegral <- function(G,B) c(0,cumsum(tail(G,-1)*diff(B)))
Ir <- rightintegral(Bvec,Bvec)

## ... and Stratonovich approximation 
Is <- (I+Ir)/2

## Plots
plot(tvec,I,type="l",col=1,ylim=range(c(I,Ir)))
lines(tvec,0.5*(Bvec^2-tvec),type="l",col=2)
lines(tvec,Ir,type="l",col=3)
lines(tvec,Is,type="l",col=4)
lines(tvec,0.5*Bvec^2,type="l",col=5)

legend(x="bottomleft",leg=c("Left","Analytical","Right","Stratonovich","Naive analytical"),
       lty="solid",col=1:5)


## -----------------------------------------------------------------------------
## Parameters and model equations
lambda <- 0.5
xi <- 2
gamma <- 1

f <- function(x) lambda*(xi-x)
g <- function(x) gamma*sqrt(abs(x))

## Initial condition
x0 <- xi

## Time points
tvec <- seq(0,100,0.01)
dt <- diff(tvec)

## Brownian motion
B <- rBM(tvec)
dB <- diff(B)

## Setup array for solution 
X <- numeric(length(tvec))
X[1] <- x0


## -----------------------------------------------------------------------------
## Simulation using the Euler stepping 
for(i in 1:(length(tvec)-1))
    X[i+1] <- X[i] + f(X[i])*dt[i] + g(X[i])*dB[i]


## -----------------------------------------------------------------------------
## Repeat with coarser time steps
t2 <- tvec[seq(1,length(tvec),10)]
dt2 <- diff(t2)

B2 <- B[seq(1,length(B),10)]
dB2 <- diff(B2)

X2 <- numeric(length(t2))
X2[1] <- x0

for(i in 1:(length(t2)-1))
    X2[i+1] <- X2[i] + f(X2[i])*dt2[i] + g(X2[i])*dB2[i]


## -----------------------------------------------------------------------------
## Plot the sample path
plot(tvec,X,type="l")

## Add line for the expectation
abline(h=xi)

lines(t2,X2,type="l",col="blue")


## -----------------------------------------------------------------------------
## Compare with the integral version
Xint <- x0+itointegral(f(X),tvec) + itointegral(g(X),B)
plot(tvec,Xint,type="l")


## -----------------------------------------------------------------------------
print(max(abs(Xint-X)))

