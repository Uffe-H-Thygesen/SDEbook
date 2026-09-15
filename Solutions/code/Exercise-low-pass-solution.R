require(SDEtools)

T <- 10
dt <- 1e-4

tv <- seq(0,T,dt)

B <- rBM(tv)

lambda <- 100

X <- euler(function(x)-lambda*x,function(x)lambda,tv,0,B)$X

Y <- stochint(X,tv)

par(mfrow=c(3,1))

plot(tv,B,type="l")
lines(tv,Y,col="red")

plot(tv,stochint(B,B),type="l")
lines(tv,stochint(Y,Y),col=2)
lines(tv,stochint(B,B,rule="c"),col=3)
lines(tv,stochint(Y,Y,rule="c"),col=4)


plot(tv,stochint(Y,B),type="l",col=1)
lines(tv,stochint(B,Y),col=2)
lines(tv,stochint(B,Y,rule="c"),col=3)
lines(tv,stochint(Y,B,rule="c"),col=4)


