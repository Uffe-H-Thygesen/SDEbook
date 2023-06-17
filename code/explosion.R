## Figure demonstrating explosion in an SDE

pdf(file="explosion.pdf",width=6,height=4)
set.seed(12345)

## Number of sample paths
N <- 3

## Time discretization
T <- 10
dt <- 0.001
tvec <- seq(0,T,dt)

nt <- length(tvec)-1

## Simulate Brownian motion
## Note: This could also have been done with SDEtools::rBM or SDEtools::rvBM
dB <- array(rnorm(N*nt,sd=sqrt(dt)),c(nt,N))
B <- apply(dB,2,function(x)c(0,cumsum(x)))

## The exploding process 
X <- tan(B)

## Threshold for explosion
Xmax <- 20

## Range for plot
Xlim <- 10

## Find time of explosion
tau <- apply(X,2,function(x)sum(cumprod(abs(x)<Xmax)))

## Set up plot
plot(c(0,max(tau))*dt,c(-1,1)*Xlim,type="n",xlab="t",ylab=expression(X[t]))

## Add sample paths
for(i in 1:N) lines(tvec[1:tau[i]],X[1:tau[i],i],lwd=2)

dev.off()
