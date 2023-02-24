## -----------------------------------------------------------------------------
## Q1


f <- function(x) r*x*(1-(x/K)^p)
g <- function(x) sigma * x

par(mfrow=c(1,2))

## arbitrary parameters as numerical examples
r <- 1
sigma <- 1

for(K in c(1,100))
    {
        ps <- c(2,1,0.25)

        p <- ps[1]
        plot(f,from=0,to=K,main=paste("K=",K))
        for(p in ps)
            plot(f,from=0,to=K,add=TRUE)
    }


## -----------------------------------------------------------------------------
## Verify the zero solution
print(f(0))
print(g(0))


## -----------------------------------------------------------------------------
r <- 0.25
sigma <- 1
K <- 1
p <- 1.5

require(SDEtools)

tv <- seq(0,100,0.01)
B <- rBM(tv)

solStab <- euler(f,g,tv,0.01,B,abs)
r <- 1
solUstab <- euler(f,g,tv,0.01,B,abs)

plot(solUstab$t,solUstab$X,type="l")
lines(solStab$t,solStab$X,col="red")


## -----------------------------------------------------------------------------
fp <-  function(x) r*(1-(p+1)*(x/K)^p)
gp <-  function(x) sigma
fXS <- function(xs) c(f(xs[1]),fp(xs[1])*xs[2])
gXS <- function(xs) c(g(xs[1]),gp(xs[1])*xs[2])

xs0 <- c(0.01*K,1)

sim <- euler(fXS,gXS,tv,xs0,B,abs)

par(mfrow=c(1,2))
plot(tv,sim$X[,1],type="l")
plot(tv,sim$X[,2],type="l",log="y")


## -----------------------------------------------------------------------------
lambdat <-  log(sim$X[,2])/tv
plot(tv,lambdat,type="l",ylim=quantile(lambdat,c(0.05,0.95),na.rm=TRUE))


## -----------------------------------------------------------------------------
r <- 0.25
sigma <- 1
K <- 1
p <- 1.5

tv <- seq(0,100,0.1)

sol <- sapply(1:100,function(i)euler(f,g,tv,0.01,B=NULL,abs)$X)
matplot(tv,sol,type="l",lty="solid")


## -----------------------------------------------------------------------------
plot(tv,apply(sol^2,1,mean),type="l")

