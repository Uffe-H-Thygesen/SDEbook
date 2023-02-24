## Simulation of the Brownian unicycle

require(SDEtools)

T <- 1
dt <- 1e-4

N <- 1000

xT <- array(0,c(N,3))

t <- seq(0,T,dt)

X <- array(0,c(N,3,length(t)))

for(i in 1:N)
    {
        B <- rBM(t)

        f <- function(x) 0*x
        g <- function(x) cbind( c(1,0,0),c(0,cos(x[1]),sin(x[1])))

        x0 <- numeric(3)

        sim <- euler(f,g,t,x0,B)
        X[i,,] <- t(sim$X)
        xT[i,] <- sim$X[length(t),]
    }

vX <- apply(X,c(2,3),var)


par(mfrow=c(3,1))

c <- c(1,1,0.5)
p <- c(1,1,2)
for(i in 1:3)
{
    plot(t,vX[i,],type="l")
    lines(t,c[i]*t^p[i])
}
