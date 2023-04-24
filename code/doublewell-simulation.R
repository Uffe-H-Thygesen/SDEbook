## Simulation of the double well
require(SDEtools)

sigma <- 0.5  ## Noise parameter

tv <- seq(0,500,0.01)
B <- rBM(tv)

f <- function(x) x*(1-x^2)
g <- function(x) sigma

sim <- euler(f,g,tv,0,B)

pdf(file="doublewell-simulation.pdf",width=6,height=5)
plot(tv,sim$X,xlab="Time",ylab="X",type="l")
dev.off()
