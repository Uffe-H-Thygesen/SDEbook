## Brownian motion on the strip x in R, y in [0,1]

T <- 1
dt <- 1e-5
n <- 3
N <- 10

## Effect of time  step on final position

X <- array(0,c(n,N))
Y <- X + 0.5

require(SDEtools)

B <- rvBM(seq(0,T,dt),N)

for(j in 1:n)
{
    for(i in 1:(T/dt))
    {
        X[j,] <- X[j,] + B[i+1,]-B[i,]
        Y[j,] <- Y[j,] + B[i+1,]-B[i,]
    }
    dt <- 2*dt
    

            
