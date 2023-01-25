require(SDEtools)

## Predator-prey dynamics
f <- function(x)
{
    C <- beta*prod(x)/(1+beta*x[1]/Cmax)
    return(c(r*x[1]*(1-x[1]/K) - C,epsilon*C-mu*x[2]))
}

g <- function(x) diag(sigma*x)

## Parameters
K <- 1
r <- 1
beta <- 2
Cmax <- 2
epsilon <- 0.15
mu <- 0.1

sigma <- rep(0.25,2)

x0 <- c(K,epsilon*K)

tvec <- seq(0,10,length=1+2^22)
B <- rvBM(tvec,n=2)

sim <- list(heun(f,g,tvec,x0,B,abs))

## matplot(sim$times,sim$X,type="l")

for(i in 1:6)
{
    I <- seq(1,length(tvec),2)
    tvec <- tvec[I]
    B <- B[I,]
    
    sim <- c(sim,list(heun(f,g,tvec,x0,B,abs)))
}

Xend <- sapply(sim,function(s) tail(s$X[,1],1))

matplot(sim[[3]]$times,sim[[3]]$X,type="l")
