## Code to simulate an autopilot designed with LQR control

set.seed(123456)

## System matrices for LQR problem
A <- matrix(c(-1,0,0,
              1,0,0,
              0,1,0)   ,nrow=3,byrow=TRUE)
B <- matrix(c(1,
              0,
              0), nrow=3)
G <- matrix(c(1,
              0,
              0), nrow=3)

## Average torque for simulation
bias <- 5

## LQR weights
Q <- diag(c(1,1,0.1))
R <- matrix(1)

require(SDEtools)

## Simulation control
T <- 100
dt <- 0.1

tv <- seq(0,T,dt)
BM <- rBM(tv,u=bias)

## Initial condition
x0 <- c(0,10,0)

## Core function to design and simulate the autopilot
autopilot <- function(Q) 
{
    ## Construct Hamiltonian matrix
    H <- rbind(cbind(A,-B %*% solve(R) %*% t(B)),cbind(-Q,-t(A)))

    ## Compute eigenvalue decomposition
    ev <- eigen(H)

    ## Identify the stable ones
    idx <- Re(ev$values) < 0

    ## Extract stable right eigenvectors
    V <- ev$vectors[,idx]
    
    ## Construct solution to Riccati equation
    n <- nrow(A)
    X <- V[ n + (1:n) ,] %*% solve( V[1:n,])

    ## Make sure that it is real and symmetric
    X <- Re(X)
    X <- 0.5*(X+t(X))
    
    ## Control gain
    F <- solve(R) %*% t(B) %*% X

    ## Closed loop system matrix
    Acl <- A - B %*% F

    ## Time-varying reference
    thetaref <- function(t) 10*(t> 0.5*tail(tv,1))

    ## Drift - notice that the second state is the error 
    f <- function(t,x) Acl %*% (x - c(0,thetaref(t),0))
    g <- function(x) G

    ## Simulate using heun method - the noise is additive, so this is
    ## consistent and more accurate than the Euler method 
    sim <- heun(f,g,tv,x0,B=BM)

    ## Add reference and control
    sim$thetaref <- thetaref(tv)
    sim$U <- - t (  (sim$X - outer(sim$thetaref,c(0,1,0)))   %*% t(F) )
    return(list(X=X,F=F,sim=sim))
}

## Design the autopilot and simulate it
pilot <- autopilot(Q)

## Plot results 
pdf(file="autopilot.pdf",width=8,height=8)

par(mfrow=c(2,1),mar=c(4,5,0.5,0.5))
plot(tv,pilot$sim$X[,2],xlab="",ylab=expression("Heading "*theta),type="l")
lines(tv,pilot$sim$thetaref,col="grey",lty="dashed")
plot(tv,pilot$sim$U,xlab="Time",ylab=expression("Rudder angle "*phi),type="l")

dev.off()
