## -----------------------------------------------------------------------------
require(SDEtools)
require(Matrix)
set.seed(123456)  ## Fix the seed so that I know what the plot looks like :) 

sigma <- 1
T <- 10

dt <- 0.01
tvec <- seq(0,T,dt)
x0 <- 0.1

fopt <- function(x) x*(1-x) - x^2
fch <- function(x) x*(1-x) - 0.5*x

g <- function(x) sigma*x

## Number of realizations
M <- 100

## Generate noise for all realizations 
B <- rvBM(tvec,n=M)

sol.opt <- sapply(1:M,function(i)euler(fopt,g,tvec,x0,B=B[,i],p=abs)$X)
sol.ch <-  sapply(1:M, function(i)euler(fch,g,tvec,x0,B=B[,i],p=abs)$X)

J.opt <- mean(sqrt(sol.opt^2))
J.ch <- mean(sqrt(0.5*sol.ch))


## -----------------------------------------------------------------------------
matplot(tvec,cbind(sol.opt[,2],sol.ch[,2]),type="l",xlab="t",ylab="X")


## -----------------------------------------------------------------------------
print(c(J.opt,J.ch))


## -----------------------------------------------------------------------------
require(SDEtools)

### Discretization of state space
Xmax <- 4
dx <- 0.01
xi <- seq(0,Xmax,dx)
xc <- xi[-1] - 0.5*diff(xi)

sigma <- 1

### Functions entering into the model

### Uncontrolled system:
D <- function(x) 1/2*sigma^2*x^2
dD <- function(x) sigma^2*x

f <- function(x) x*(1-x)
advection <- function(x) f(x) - dD(x)

G0 <- fvade(advection,D,xi,'r')

### Effect of the fishing: The "generator" d/dx
G1 <- fvade(function(x)-1,function(x)0,xi,'r')

ubound <- c(0,rep(100,nrow(G1)-1))

k <- function(u) sqrt(u)

vbar <- c(Inf,rep(0.01,nrow(G0)-1))
hack <- function(dV) pmax(-dV,vbar)
uopt <- function(dV) 1/4/hack(dV)^2

sol <- PolicyIterationSingular(G0,G1,k,uopt,do.minimize = FALSE)


## -----------------------------------------------------------------------------
par(mfrow=c(2,1))
plot(xc,sol$V - approx(xc,sol$V,1)$y,xlab="x",ylab="V")
lines(xc,0.5*log(xc))
      
plot(xc,sol$u,ylim=c(0,max(xc)^2),xlab="x",ylab="u")
lines(xc,xc^2)


## -----------------------------------------------------------------------------
pv <- c(0.5,1,2)

sols <- list()


for(p in pv) 
{
    f <- function(x) x*(1-x^p)
    advection <- function(x) f(x) - dD(x)
    G <- fvade(advection,D,xi,'r')

    sol <- PolicyIterationSingular(G,G1,k,uopt,do.minimize = FALSE)
    
    sols[[length(sols)+1]] <- list(p=p,sol=sol)
}


## -----------------------------------------------------------------------------
par(mfrow=c(2,1))
matplot(xc,cbind(sols[[1]]$sol$V,sols[[2]]$sol$V,sols[[3]]$sol$V),
	type="l",xlab="x",ylab="V")
matplot(xc,cbind(sols[[1]]$sol$u,sols[[2]]$sol$u,sols[[3]]$sol$u),
	type="l",xlab="x",ylab="u",ylim=c(0,20))

