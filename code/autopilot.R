set.seed(123456)

A <- matrix(c(-1,0,0,
              1,0,0,
              0,1,0)   ,nrow=3,byrow=TRUE)
B <- matrix(c(1,
              0,
              0), nrow=3)
G <- matrix(c(1,
              0,
              0), nrow=3)

bias <- 5

Q <- diag(c(1,1,0.1))

R <- matrix(1)

require(SDEtools)

T <- 100
dt <- 0.1

tv <- seq(0,T,dt)
BM <- rBM(tv,u=bias)

x0 <- c(0,10,0)


autopilot <- function(Q) 
  {
    H <- rbind(cbind(A,-B %*% solve(R) %*% t(B)),cbind(-Q,-t(A)))

    ev <- eigen(H)

    idx <- Re(ev$values) < 0

    V <- ev$vectors[,idx]
    
    n <- nrow(A)
    X <- V[ n + (1:n) ,] %*% solve( V[1:n,])
    X <- Re(X)
    X <- 0.5*(X+t(X))
    
    F <- solve(R) %*% t(B) %*% X

    evs <- ev$values[idx]

    Acl <- A - B %*% F

  thetaref <- function(t) 10*(t> 0.5*tail(tv,1))
  f <- function(t,x) Acl %*% (x - c(0,thetaref(t),0))
  g <- function(x) G

  sim <- heun(f,g,tv,x0,B=BM)

  sim$thetaref <- thetaref(tv)
  sim$U <- - t (  (sim$X - outer(sim$thetaref,c(0,1,0)))   %*% t(F) )
  return(list(X=X,F=F,sim=sim))
}



pilot <- autopilot(Q)

pdf(file="autopilot.pdf",width=8,height=8)

par(mfrow=c(2,1),mar=c(4,5,0.5,0.5))
plot(tv,pilot$sim$X[,2],xlab="",ylab=expression("Heading "*theta),type="l")
lines(tv,pilot$sim$thetaref,col="grey",lty="dashed")
plot(tv,pilot$sim$U,xlab="Time",ylab=expression("Rudder angle "*phi),type="l")

dev.off()
