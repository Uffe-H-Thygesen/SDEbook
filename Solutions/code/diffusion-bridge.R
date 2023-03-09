require(SDEtools)

sigma <- 0.5

f <- function(x) -x^3+x
g <- function(x) sigma
D <- function(x) 0.5*g(x)^2 

T <- 100
xi <- seq(-3,3,length=200)
xc <- 0.5*(xi[-1]+xi[-length(xi)])
nx <- length(xc)

G <- fvade(f,D,xi,'r')

x0 <- -1
xT <- 1

phi0 <- numeric(length(xi)-1)
phi0[sum(xi<x0)] <- 1

psiT <- numeric(length(xi)-1)
psiT[sum(xi<xT)] <- 1

nt <- 100
dt <- T/nt
tv <- seq(0,T,dt)

P <- expm(G*dt)

PSI <- array(NA,c(nx,nt+1))
PHI <- array(NA,c(nt+1,nx))
PI <- array(NA,c(nt+1,nx))

PHI[1,] <- phi0
PSI[,nt+1] <- psiT

for(i in 1:nt) PHI[i+1,] <- as.numeric(PHI[i,] %*% P)
for(i in nt:1) PSI[,i] <- as.numeric(P %*% PSI[,i+1])
for(i in 1:(nt+1)) PI[i,] <- PHI[i,] * as.numeric(PSI[,i])
for(i in 1:(nt+1)) PI[i,] <- PI[i,] / sum(PI[i,])

image(tv,xc,t(apply(PI,1,cumsum)))

EX <- PI %*% xc
EX2 <- PI %*% xc^2 
sX <- sqrt(EX2-EX^2)

lines(tv,EX)
lines(tv,EX+sX,lty="dashed")
lines(tv,EX-sX,lty="dashed")
