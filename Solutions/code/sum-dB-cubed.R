## Verify by simulation that the sum of cubed increments in BM vanish
## as the mesh becomes finer
require(SDEtools)

t <- 1
tv <- seq(0,t,length=2^20+1)

B <- rBM(tv) 

## Number of grid doublings
n <- 10

sB3 <- numeric(n)

dt <- (tv[2]-tv[1])*2^(0:(n-1))

for(i in 1:n)
{
    sB3[i] <- sum(abs(diff(B)^3))
    B <- B[seq(1,length(B),2)] 
}

plot(dt,sB3,log="xy",xlab="Time step",ylab="sum dB^3")

lines(dt,2*sqrt(dt),lty="dashed")
