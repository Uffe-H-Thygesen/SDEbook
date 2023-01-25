require(SDEtools)

lambda <- 1000
sigma <- 1

T <- 10
dt <- 1e-5

tv <- seq(0,T,dt)
B <- rBM(tv)

fC <- function(x)c(-lambda*x[1],-sin(x[2]) + sigma*cos(x[2])*x[1])
gC <- function(x)c(lambda,0)

x0 <- c(0,0)

simCe <- euler(fC,gC,tv,x0,B)
simCh<- heun(fC,gC,tv,x0,B)

plot(simCh$t,simCh$X[,2],type="l")
lines(simCe$t,simCe$X[,2],type="l",col="red")

f <- function(x) -sin(x) 
g <- function(x) cos(x)

simI <- euler(f,g,tv,0,B)
simS <- heun(f,g,tv,0,B)

lines(simI$t,simI$X,type="l",col="blue")
lines(simS$t,simS$X,type="l",col="green")

err <- cbind(simCe$X[,2],simI$X,simS$X) - simCh$X[,2]

var(err)

