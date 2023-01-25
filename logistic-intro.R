## Simulation of logistic growth model

set.seed(12345)

r <- 1
K <- 1
sigma <- 0.1

f <- function(x) r*x*(1-x/K)
g <- function(x) sigma*x

x0 <- 0.01 * K

T <- 20/r
dt <- 0.001/r

tv <- seq(0,T,dt)

require(SDEtools)

sim0 <- euler(f,function(x)0,tv,x0,p=function(x)abs(x))
sim1 <- euler(f,g,tv,x0,p=function(x)abs(x))
sim2 <- euler(f,g,tv,x0,p=function(x)abs(x))

ylim <- c(0,max(c(sim0$X,sim1$X,sim2$X)))

pdf(file="logistic-intro.pdf",height=3,width=6)
par(mar=c(4.5,4,0.2,0.2))
plot(sim0$t,sim0$X,type="l",ylim=ylim,xlab="t",ylab="X")
lines(sim1$t,sim1$X)
lines(sim2$t,sim2$X,lty="dotted")
dev.off()
