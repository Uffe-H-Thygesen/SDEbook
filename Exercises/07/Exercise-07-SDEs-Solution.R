## -----------------------------------------------------------------------------
require(SDEtools)

## Define drift and noise terms
f <- function(x) c(x[2],sin(x[1]) - lambda*x[2])
g <- function(x) c(0,sigma)

## System parameters
lambda <- 0.1
sigma <- 0.1

## Simulation parameters
dt <- 0.01
x0 <- c(0,0)
T <- 1000

tv <- seq(0,T,dt)

## Perform the simulation with the Euler method
B <- rBM(tv)
sol <- euler(f,g,tv,x0,B)

## Plot also the last 10 percent
I <- tv>0.9*T

par(mfrow=c(2,2))
plot(sol$times,sol$X[,1],type="l",xlab="t",ylab="x")
hist(sol$X[,1],main="",xlab="x",freq=FALSE)
plot(sol$times[I],sol$X[I,1],type="l",xlab="t",ylab="x")
hist(sol$X[I,1],main="",xlab="x",freq=FALSE)


## -----------------------------------------------------------------------------
sigma <- 1/2
sol <- euler(f,g,tv,x0,B)
par(mfrow=c(2,2))
plot(sol$times,sol$X[,1],type="l",xlab="t",ylab="x")
hist(sol$X[,1],main="",xlab="x",freq=FALSE)
plot(sol$times[I],sol$X[I,1],type="l",xlab="t",ylab="x")
hist(sol$X[I,1],main="",xlab="x",freq=FALSE)


## -----------------------------------------------------------------------------
lambda <- mu <- 1

sim <- function(sigma,tv,x0=mu,B=NULL)
    {
        g <- function(x) sigma*sqrt(abs(x))
        fI <- function(x) lambda*(mu-x)
        fS <- function(x) fI(x) - sigma^2/4

        if(is.null(B)) B <- rBM(tv)

        simI <- euler(f=fI,g=g,times=tv,x0=x0,B=B,p=abs)
        simS <- heun (f=fS,g=g,times=tv,x0=x0,B=B,p=abs)
        
        return(cbind(XI = simI$X,XS = simS$X))
    }


## -----------------------------------------------------------------------------
sigmas <- seq(0.25,1.75,length=4)
tv <- seq(0,100,0.01)
B <- rBM(tv)
sols <- lapply(sigmas,function(s)sim(s,tv,B=B))


## -----------------------------------------------------------------------------
par(mfrow=c(4,2))
xmax <- 4
xi <- seq(0,xmax,0.1)
for(sol in sols)
    {
        plot(tv,sol[,1],type="l",xlab="t",ylab="x")
        lines(tv,sol[,2],col="red")
        hist(pmin(sol[,2],xmax),breaks=xi,freq=FALSE,xlab="x",main="")
    }


## -----------------------------------------------------------------------------
sigma <- 1.25
tv1 <- seq(0,100,0.01)
B1 <- rBM(tv1)
sol1 <- sim(sigma,tv1,B=B1)

I <- seq(1,length(B1),10)
B2 <- B1[I]
tv2 <- tv1[I]

sol2 <- sim(sigma,tv2,B=B2)

## Plot only the tail, for increased clarity
matplot(tv1,sol1,type="l",xlim=c(90,100))
matplot(tv2,sol2,type="l",add=TRUE)


## -----------------------------------------------------------------------------
## Difference between Ito and Stratonovich on the fine grid 
print(mean(abs(sol1[,1]-sol1[,2])))

## Difference between Ito and Stratonovich on the coarse grid 
print(mean(abs(sol2[,1]-sol2[,2])))

## Difference between fine and coarse grid for Ito 
print(mean(abs(sol1[I,1]-sol2[,1])))

## Difference between fine and coarse grid for Ito 
print(mean(abs(sol1[I,2]-sol2[,2])))

