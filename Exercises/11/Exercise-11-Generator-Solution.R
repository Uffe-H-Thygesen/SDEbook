## ----message=FALSE------------------------------------------------------------
require(SDEtools)
require(bvpSolve)
require(MASS)

f <- function(x) -x
g <- function(x) sqrt(2)

x0 <- 0.5
l <- 2

## "Reasonable" time step 
dt <- 1e-3


## -----------------------------------------------------------------------------
## Simulate with the Euler method until exit from domain
sim <- function(dummy)
{
    x <- x0
    t <- 0
    r <- 0
    while(abs(x)<l)
    {
        t <- t + dt
        x <- x + f(x)*dt + g(x)* rnorm(1,sd=sqrt(dt))
        r <- r + x^2 *dt 
    }

    return(c(t=t,x=x,r=r))
}


## -----------------------------------------------------------------------------
## Obtain N samples
N <- 1000
res <- sapply(1:N,sim)

print(paste("Mean time to exit: ",mean(res[1,])))
print(paste("Probability of exit to the right: ",mean(res[2,]>0)))


## -----------------------------------------------------------------------------
## Solving the BVP with bvpSolve
fun <- function(x,y,p)
    list(c(y[2],-2/g(x)/g(x)*f(x)*y[2]))

x <- seq(-l,l,length=201)
sol <- bvpshoot(yini = c(0,NA),yend=c(1,NA),x=x,func=fun,guess=0)


## -----------------------------------------------------------------------------
## Compute the solution from the scale function, using numerical integration
phi <- function(x) exp(integrate(function(y)-2*f(y)/g(y)/g(y),lower=0,upper=x)$value)
phiv <- function(x) sapply(x,phi)  # Vectorize
s <- function(x) integrate(phiv,lower=0,upper=x)$value
sv <- function(x) sapply(x,s)      # Vectorize

sm <- s(-l)  # Boundary points 
sp <- s(+l)

ss <- function(x) (sv(x)-sm)/(sp-sm)   # Scaled scale function (sic)


## -----------------------------------------------------------------------------
## Discretize the generator and find its null space
xi <- seq(-l,l,length=201)
xc <- xi[-1] - 0.5*diff(xi)
G <- fvade(f,function(x)0.5*g(x)*g(x),xi,'e')

sG <- Null(t(G))

coefs <- solve(sG[c(1,nrow(sG)),],c(0,1))

h <- sG %*% coefs


## -----------------------------------------------------------------------------
## bvpsolve
plot(sol[,1],sol[,2],xlab="x",ylab="Probability of exit right")

## Monte carlo
points(x0,mean(res[2,]>0))
lines(c(x0,x0),mean(res[2,]>0)+c(-2,2)*sqrt(var(res[2,]>0)/ncol(res)))

## Scale function by numerical integration 
plot(ss,from=-l,to=l,add=TRUE,col="red")

## With the generator
lines(xc,h[2:(length(h)-1)],col="green")


## -----------------------------------------------------------------------------
## Expected time to exit
fun <- function(x,y,p)
    list(c(y[2],-2/g(x)/g(x)*(f(x)*y[2]+1)))

x <- seq(-l,l,length=201)
sol.bvp <- bvpshoot(yini = c(0,NA),yend=c(0,NA),x=x,func=fun,guess=0)


## -----------------------------------------------------------------------------
G <- fvade(f,function(x)0.5*g(x)*g(x),xi,'a')
k.fvade <- solve(G,rep(-1,nrow(G)))


## -----------------------------------------------------------------------------
## Plot bvp solution 
plot(sol.bvp[,1],sol.bvp[,2],type="l",xlab="x",ylab="Expected time to exit")

## Add solution from generator
lines(xc,k.fvade,col="green")

## Add Monte Carlo; add confidence interval 
points(x0,mean(res[1,]))
lines(rep(x0,2),mean(res[1,]) + c(2,-2)*sqrt(var(res[1,])/N))


## -----------------------------------------------------------------------------
## Cumulated reward (x squared)
fun <- function(x,y,p)
    list(c(y[2],-2/g(x)/g(x)*(f(x)*y[2]+x^2)))

x <- seq(-1,1,0.01)
sol.bvp <- bvpshoot(yini = c(0,NA),yend=c(0,NA),x=x,func=fun,guess=0)


## Compute solution from generator
G <- fvade(f,function(x)0.5*g(x)*g(x),xi,'a')
k.fvade <- solve(G,-xc^2)


## -----------------------------------------------------------------------------
plot(sol.bvp[,1],sol.bvp[,2],type="l",xlab="x",ylab="Expected cumulated reward")
lines(xc,k.fvade,col="green")

## Add Monte Carlo; add confidence interval 
points(x0,mean(res[3,]))
lines(rep(x0,2),mean(res[3,]) + c(2,-2)*sqrt(var(res[3,])/N))

