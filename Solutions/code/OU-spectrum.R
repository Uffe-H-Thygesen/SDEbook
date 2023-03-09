## Compute the slowest modes (forward and backward) for the Ornstein-Uhlenbeck process

require(SDEtools)
require(RSpectra)

f <- function(x) -x
D <- function(x) 1

xi <- seq(-4,4,0.01)
xc <- 0.5*(head(xi,-1)+tail(xi,-1))

G <- fvade(f,D,xi,'r')


evsBWD <- eigs(G,k=4,sigma=1e-8)
evsFWD <- eigs(t(G),k=4,sigma=1e-8)

evsBWD$values <- Re(evsBWD$values)
evsFWD$values <- Re(evsFWD$values)
evsBWD$vectors <- Re(evsBWD$vectors)
evsFWD$vectors <- Re(evsFWD$vectors)


## Coefficients in the Hermite polynomials
H <- matrix(c(0,0, 0, 1,
              0,0, 1, 0,
              0,1, 0,-1,
              1,0,-3, 0),nrow=4,byrow=TRUE)

X <- outer(3:0,xc,function(p,x)x^p)
hp <- H %*% X

rho <- dnorm(xc)

par(mfrow=c(4,2))
for(i in 4:1)
{
    j <- 5-i
    
    ev <- hp[j,]
    lambdatheo <- -j+1
    lambdanum  <- evsBWD$values[i]
    
    plot(xc,ev,type="l",main=paste("BWD. Theo:",lambdatheo,"Num: ",lambdanum))
    k <- sum(evsBWD$vectors[,i] * ev)/sum(ev^2)
    lines(xc,evsBWD$vectors[,i] / k,col="red")
    
    ev <- rho * hp[j,]
    plot(xc,ev,type="l",main="FWD")
    k <- sum(evsFWD$vectors[,i] * ev)/sum(ev^2)
    lines(xc,evsFWD$vectors[,i] / k,col="red")
}
