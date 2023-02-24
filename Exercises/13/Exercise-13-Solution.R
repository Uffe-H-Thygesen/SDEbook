require(fields)
require(SDEtools)
require(MASS)
require(expm)

rm(list=ls())
graphics.off()

## Simulate the optimal and the "constant harvest rate" control strategy
sigma <- 1

ustar <- function(x) x^2
fopt <- function(x) x*(1-x) - ustar(x)

uch <- function(x) 0.5*x
fch <- function(x) x*(1-x) - uch(x)

g <- function(x) sigma*x

T <- 10
dt <- 0.01
tvec <- seq(0,T,dt)
x0 <- 1
alpha <- 0.0 # Cost of harvest rate

B <- rBM(tvec)
solopt <- euler(fopt,g,tvec,x0,B)
solch <- euler(fch,g,tvec,x0,B)

plot(solopt$times,solopt$X,type="l",col="blue",xlab="Time",ylab="X")
lines(solch$times,solch$X, type="l",col="red")

dev.copy2pdf(file="simFish.pdf")

## Discretize the HJB equation 
Xmax <- 3
dx <- 0.03
xi <- seq(0,Xmax,dx)
xc <- xi[-1] - 0.5*diff(xi)

Umax <- 1

f0 <- function(x) x*(1-x)
b <- function(x) -1

D <-  function(x) 0.5*sigma^2*x^2
Dp <- function(x) sigma^2*x

L0 <- fvade(function(x) f0(x)-Dp(x),D,xi,'r')
L1 <- fvade(b,function(x)0,xi,'r')

ddx <- fvade(function(x)1,function(x)0,xi,'a')


T <- 10
dt <- 0.001

tvec <- seq(0,T,dt)

V <- array(NA,c(length(xc),length(tvec)))
U <- V

V[,length(tvec)] <- 0
U[,length(tvec)] <- 0

## Time-step the unctrolled part using expm rather than Euler. This
## allows larger time steps
P <- expm(L0*dt)

for(i in (length(tvec)-1):1)
    {
        ## Compute optimal control at each point in state space
        dVdx <- - L1 %*% V[,i+1]
        
        ## Optimal control without penalty on harvest rate
        ## Truncate at Umax 
        U[,i] <- pmin(1/4/(dVdx^2),4*xc)

        ## Optimal control with penalty on harvest rate (CHECK THIS)
        ## U[,i] <- pmin(xc/(2*dVdx*xc+2*alpha)^2,Umax)

        ## Update HJB equation
        V[,i] <- P %*% V[,i+1] + U[,i] * L1 %*% V[,i+1] * dt + sqrt(U[,i]) * dt - alpha*U[,i]*dt
    }

X11()
par(mfrow=c(2,1))

## Supsample time to get reasonably sized plots
I <- seq(1,length(tvec),10)
image.plot(tvec[I],xc,t(V[,I]),main="V",xlab="Time",ylab="x")
image.plot(tvec[I],xc,t(U[,I]),main="U",xlab="Time",ylab="x")

dev.copy2pdf(file="OptimalControl.pdf")


## Extract pieces from the solution
X11()
par(mfrow=c(2,1))
plot(tvec,V[round(length(xc)/2),],type="l",xlab="Time",ylab="V(0.5,t)")
plot(xc,U[,1],xlab="x",ylab="Optimal control",type="l")
dev.copy2pdf(file="slice.pdf")

