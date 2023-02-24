require(SDEtools)

r <- 1

x0 <- 0

sigma <- 0.5

tv <- seq(0,500,0.01)
B <- rBM(tv)

f <- function(x) r*x*(1-x^2)
g <- function(x) sigma

sim <- euler(f,g,tv,x0,B)

pdf(file="doublewell-simulation.pdf",width=6,height=5)
plot(tv,sim$X,xlab="Time",ylab="X",type="l")
dev.off()
