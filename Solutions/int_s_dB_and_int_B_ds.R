## Monte Carlo simulation to verify the variance and covariance of int s dB and int B ds

## Pick an arbitrary terminal time, and a reasonably fine resolution
t <- 5
tv <- seq(0,t,length=1001)

require(SDEtools)

## Helper function to compute one realization of the two integrals
fun <- function(i)
{
    Bv <- rBM(tv)
    I1 <- stochint(tv,Bv)[length(tv)]
    I2 <- stochint(Bv,tv)[length(tv)]

    return(c(I1,I2))
}

## Generate 10000 samples (maybe a bit overkill)
Is <- sapply(1:10000,fun)

## Display variance/covariance
var(t(Is))

## Compare with the theoretical expectation
array(c(1/3,1/6,1/6,1/3),c(2,2))*t^3
