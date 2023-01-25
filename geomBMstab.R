T <- 10
dt <- 0.01

tvec <- seq(0,T,dt)

nt <- length(tvec)-1

dB <- rnorm(nt,mean=0,sd=sqrt(dt))

B <- c(0,cumsum(dB))

r <- 0.25
s <- 1

s2 <- s^2

X <- exp( (r-0.5*s2)*tvec + s* B)

plot(tvec,X,type="l")

lines(tvec,exp(r*tvec))

lines(tvec,exp(qnorm(
