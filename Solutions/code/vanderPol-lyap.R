## Assess the Lyapunov exponent of the van der Pol oscillator

require(SDEtools)
## set.seed(123456)

## Define the system
fx <- function(x,v) v
fv <- function(x,v) v*(mu - x^2)-x
gx <- function(x,v) 0
gv <- function(x,v) sigma

f <- function(xv) c(fx(xv[1],xv[2]),fv(xv[1],xv[2]))
g <- function(xv) c(0,gv(xv[1],xv[2]))

## Jacobians for sensitivity equations
jacf <- function(xv) matrix(c(0,1,-2*mu*xv[1]*xv[2]-1,mu*(1-xv[1]^2)),nrow=2,byrow=TRUE)
jacg <- function(xv) array(0,c(2,2))

## Total dynamics of X,V, and S
ff <- function(xvs) c(f(xvs[1:2]),as.numeric(jacf(xvs[1:2]) %*% array(xvs[3:6],c(2,2))))
gg <- function(xvs) c(g(xvs[1:2]),as.numeric(jacg(xvs[1:2]) %*% array(xvs[3:6],c(2,2))))

## Parameters
sigma <- 0.5
mu <- 1
T <- 100
dt <- 0.01

times <- seq(0,T,dt)
B <- rBM(times)

## Initial conditions
xv0 <- c(0,0)
xvs0 <- c(xv0,c(1,0,0,1))

sim <- heun(ff,gg,times,xvs0,B)

## Check sensitivity equations by comparison with perturbed initial state
epsilon <- 1e-8
dxv0 <- rnorm(2)*epsilon
sim1 <- heun(ff,gg,times,xvs0+c(dxv0,0,0,0,0),B)

graphics.off()
par(mfrow=c(2,1))

## Plot trajectories - the original one, and the one with the perturbed initial state
matplot(times,sim$X[,1:2],type="l",lty=1,xlab="Time",ylab="X,V")
matplot(times,sim1$X[,1:2],type="l",lty=2,add=TRUE,col=3:4)

## Plot the difference between the two trajectories, and the what the sensitivity predicts
matplot(times,sim1$X[,1:2]-sim$X[,1:2],type="l",lty=1,xlab="Time",ylab="dX,dV")
matplot(times,t(apply(sim$X[,3:6],1,function(s) as.numeric(array(s,c(2,2)) %*% dxv0))),
        add=TRUE,type="l",lty=2,col=3:4)

dev.new()

dx <- apply(sim$X[,1:2] - sim1$X[,1:2],1,function(x)sqrt(sum(x^2)))
svs <- apply(sim$X[,3:6],1,function(s) base::norm(array(s,c(2,2)),"2"))

plot(times,log(svs)/times,type="l",ylim=c(-0.1,0.1),xlab="Time",ylab=expression(lambda[t]))
lines(times,log(dx/epsilon)/times,lty=3)
