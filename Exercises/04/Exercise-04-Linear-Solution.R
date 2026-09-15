## -----------------------------------------------------------------------------
require(SDEtools)

## System parameters
m <- 1    # [kg]
k <- 0.5  # [N/m]
c <- 0.2  # [N*s/m]
sigma <- 10 # [N sqrt(s)]

## Q 1.1
A <- array(c(0,-k/m,1,-c/m),c(2,2))
G <- array(c(0,sigma/m),c(2,1))


## Simulation parameters
T <- 1000   
dt <- 0.01

## Setup arrays
tvec <- seq(0,T,dt)
P <- array(0,c(2,2,length(tvec)))  # This also computes the solution of the Lyapunov equation 
X <- array(0,c(2,length(tvec)))

## Simulate sample path of Brownian motion 
B <- rBM(tvec)
dB <- diff(B)

## Main time loop, Euler stepping the SDE and the Lyapunov equation 
for(i in 1:(length(tvec)-1))
{
    X[,i+1] <- X[,i] + A %*% X[,i] * dt + G * dB[i]
    P[,,i+1] <- P[,,i] + (A %*% P[,,i] + P[,,i] %*% t(A) + G %*% t(G)) * dt
}

## Q 2.2: Plot the sample path; first part only for clarity
par(mfrow=c(2,1))
plot(tvec,X[1,],type="l",xlim=c(0,100),xlab="Time",ylab="Position")
plot(tvec,X[2,],type="l",xlim=c(0,100),xlab="Time",ylab="Velocity")


## -----------------------------------------------------------------------------

Pinf <- P[,,length(tvec)]
Pinf2 <- lyap(A,G%*%t(G))

## Print the empirical covariance and compare with the two solutions of the Lyapunov equation 
print(cov(t(X)))
print(Pinf)
print(Pinf2)


## -----------------------------------------------------------------------------
print(Ekin <- 0.5*m*var(X[2,]))
print(Epot <- 0.5*k*var(X[1,]))

## ... add the analytical predictions:

print(0.5*m*Pinf2[2,2])
print(0.5*k*Pinf2[1,1])


## -----------------------------------------------------------------------------
acf(X[1,],lag.max=50/dt)

ivec <- seq(0,50/dt,1)
tvec <- ivec * dt

require(Matrix)
rhovec <- sapply(tvec,function(t) (Pinf2 %*% expm(t(A)*t) )[1,1])
lines(ivec,rhovec/rhovec[1],col="red",lwd=3)


## -----------------------------------------------------------------------------
I <- diag(c(1,1))
H <- function(w) (solve(1i*w*I-A) %*% G)[1]
omegaR <- abs(Im(eigen(A)$values[1]))
ws <- omegaR*10^(seq(-0.5,1,length=101))
Hs <- sapply(ws,H)
par(mfrow=c(3,1))
plot(ws,abs(Hs),log="xy",type="l",ylab="|H|")
plot(ws,Arg(Hs),log="x",type="l",ylab="|H|")
plot(ws,abs(Hs)^2,type="l",log="xy",ylab="Var.spec.")
lines(ws,rep(sigma^2,length(ws)),lty="dashed")


## -----------------------------------------------------------------------------
## Specific examples of system parameters 
lambda <- 1
sigma <- 1
x <- 1

EXt <- function(t) x*exp(-lambda*t)

## Differential lyapunov equation
Lyap <- function(V) -2*lambda*V + sigma^2

## Analytical solution of the equation
VXt <- function(t) sigma^2/2/lambda*(1-exp(-2*lambda*t))

## Limit as time goes to infinity, assuming a<0
VXinf <- sigma^2/2/lambda

## Test that it is an equilibrium point for the Lyapunov equation
print(Lyap(VXinf))

plot(EXt,from=0,to=10,lty="dashed",ylim=c(-1,1)*1.5)
plot(function(t) EXt(t)+sqrt(VXt(t)),from=0,to=10,add=TRUE,lty="dotted")
plot(function(t) EXt(t)-sqrt(VXt(t)),from=0,to=10,add=TRUE,lty="dotted")

abline(h= sqrt(VXinf),lty="dotdash")
abline(h=-sqrt(VXinf),lty="dotdash")

times <- seq(0,10,0.01)
sim <- euler(f=function(x)-lambda*x,g=function(x) sigma,times=times,x0=x)
lines(sim$times,sim$X)

