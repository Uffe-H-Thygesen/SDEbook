## Simulation of a wind turbine

require(SDEtools)

set.seed(123)  ## For reproucible results

## Initial condition 
x0 <- c(0.5,0,0)   ## State variables: Force, position, velocity

## system parameters 
sigma <- 1         ## Noise on the force
lambda <- 0.5      ## Relaxation rate of force
k <- 1             ## Spring constant
mu <- 0.5          ## Damping
fbar <- 3          ## Mean force

## Construct system matrices
A <- matrix(c(-lambda,0,  0,
              0,      0,  1,
              1,     -k,-mu),nrow=3,byrow=TRUE)

G <- c(sigma,0,0)

xbar <- c(fbar,0,0)

## Define system functions in dX = f(X)*dt + g(X)*dB
f <- function(fxv) A %*% ( fxv - xbar) 
g <- function(fxv) G

## Time grid
tv <- seq(0,30,0.01)

## Simulate with and without noise 
sim <- euler(f,g,tv,x0)
sim0 <- euler(f,g,tv,x0,B=rep(0,length(tv)))

## Plots
pdf("wind.pdf",width=4,height=4)
par(mfrow=c(3,1),mar=c(2,5,2,2))

plot(sim$t,sim$X[,1],type="l",ylim=c(0,max(sim$X[,1])),ylab="Force",xlab="Time")
lines(sim0$t,sim0$X[,1],lty="dashed")
par(mar=c(3,5,1,2))

plot(sim$t,sim$X[,2],type="l",ylab="Position",xlab="Time",
     ylim=c(-2,2))
lines(sim0$t,sim0$X[,2],lty="dashed")

par(mar=c(4,5,0,2))
plot(sim$t,sim$X[,3],type="l",ylab="Velocity",xlab="Time")
lines(sim0$t,sim0$X[,3],lty="dashed")

dev.off()

## Extend simulation to estimate variance
tv <- seq(0,1000,0.01)
sim <- euler(f,g,tv,x0)

burnin <- 10

Pemp <- var(sim$X[tv>=burnin,])    ## Empirical variance
Ptheo <- lyap(A,outer(G,G))        ## Theoretical prediction
