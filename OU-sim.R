## Simulation of the Ornstein-Uhlenbeck process
require(SDEtools)
set.seed(1234567)

lambda <- 1
sigma <- 1
x0 <- 1

T <- 10
dt <- 0.01

tv <- seq(0,T,dt)

f <- function(x)-lambda*x
g <- function(x)sigma

sim <- euler(f,g,tv,x0)

mx <- x0*exp(-lambda*tv)
sx <- sqrt(sigma^2/2/lambda*(1-exp(-2*lambda*tv)))

pdf(file="OU-sim.pdf",width=6,height=3.5)
par(mar=c(4.5,4,0.1,0.1))
plot(tv,sim$X,type="l",xlab="t",ylab=expression(X[t]))

lines(tv,mx,lty="dashed")
lines(tv,mx+sx,,lty="dotted")
lines(tv,mx-sx,lty="dotted")
dev.off()
