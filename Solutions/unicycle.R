## Simulation of the Brownian unicycle

require(SDEtools)

N <- 1000
T <- 0.001
dt <- 1e-5

times <- seq(0,T,dt)

f <- function(z) numeric(3)
g <- function(z) array(c(1,0,0,0,cos(z[1]),sin(z[1])),c(3,2))

x0 <- numeric(3)
    
sim <- function()
{
    Z <- euler(f,g,times,x0)
    return(Z$X)
}

Z <- array(NA,c(length(times),3,N))

for(i in 1:N) Z[,,i] <- sim()

plot(Z[length(times),2,],Z[length(times),3,])
