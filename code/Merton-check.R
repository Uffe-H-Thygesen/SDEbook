require(SDEtools)

r <- 0.
alpha <- 0.1
sigma <- 0.4
gamma <- 0.75

x <- 1
w <- 1

p <- (alpha-r)/gamma/sigma^2

T <- 1000
dt <- 0.01

times <- seq(0,T,dt)

B <- rBM(times)

X <- x*exp((alpha-sigma^2/2)*times + sigma*B)
W <- w*exp((r+p*(alpha-r)-sigma^2*p^2/2)*times + sigma*p*B)

N <- p*W/X
Nc <- p*w/x +
    (p-1)*(alpha-r)*(gamma+1)*itointegral(N,times) +
    sigma*(p-1)*itointegral(N,B)

par(mfrow=c(2,2))
plot(times,X,type="l")
plot(times,W,type="l")
plot(times,N,type="l")
lines(times,Nc,col=2)




