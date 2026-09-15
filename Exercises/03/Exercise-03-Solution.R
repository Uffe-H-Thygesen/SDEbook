dt <- 2^-20
t <- seq(0,1,dt)
dB <- rnorm(length(t)-1,sd=sqrt(dt))
B <- c(0,cumsum(dB))
plot(t,B,type="l")
dtvec <- QV <- TV <- numeric(10)

for(i in 1:10)
{
    dtvec[i] <- dt
    TV[i] <- sum(abs(dB))
    QV[i] <- sum(dB^2)
    B <- B[seq(1,length(B),2)]
    dB <- diff(B)
    dt <- 2* dt
}


par(mfrow=c(1,2))

plot(dtvec,TV,log="xy")
lines(dtvec,sqrt(2/pi/dtvec))

plot(dtvec,QV,log="xy")
abline(h=1)


N <- 1000
dt <- 2^-10
t <- seq(0,1,dt)

S <- numeric(N)

for(i in 1:N)
{
    dB <- rnorm(length(t)-1,sd=sqrt(dt))
    B <- c(0,cumsum(dB))
    M <- exp(B-0.5*t)
    S[i] <- max(M)
    
}
