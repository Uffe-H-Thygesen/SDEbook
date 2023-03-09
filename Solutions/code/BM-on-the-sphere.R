## Verification of the expectation of BM on the sphere

require(SDEtools)

n <- 4
N <- 100
T <- 1
dt <- 0.01

times <- seq(0,T,dt)

I <- diag(n)

fS <- function(x) 0*x
fI <- function(x) (1-n)/2*x
g <- function(x) I - outer(x,x)/sum(x*x)

x0 <- c(1,rep(0,n-1))

solI <- array(NA,c(length(times),n,N))
solS <- array(NA,c(length(times),n,N))

for(i in 1:N)
{
    B <- rvBM(times,n=n)
    
    solI[,,i] <- euler(fI,g,times,x0,B)$X
    solS[,,i] <- heun(fS,g,times,x0,B)$X
}

dISm <- apply(solI-solS,1,mean)
dISv <- apply(apply(solI-solS,c(1,2),var),1,mean)

ucl <- dISm + sqrt(dISv)
lcl <- dISm - sqrt(dISv)

graphics.off()

par(mfrow=c(2,2))

plot(times,dISm,ylim=range(c(ucl,lcl)))
lines(times,ucl)
lines(times,lcl)

norm2I <- apply(solI,c(1,3),function(x)sum(x^2))
norm2S <- apply(solS,c(1,3),function(x)sum(x^2))

matplot(times,norm2I,type="l",col=1)
matplot(times,norm2S,type="l",col=2,add=TRUE)

EX1I <- apply(solI[,1,],1,mean)
EX1S <- apply(solS[,1,],1,mean)
matplot(times,EX1I,type="l",col=1)
matplot(times,EX1S,type="l",col=2,add=TRUE)

lines(times,exp((1-n)/2*times),col=3)
