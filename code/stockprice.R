# Figure generating stock prices
# For introduction


Tstop <- 10
dt <- 0.01

N <- round(Tstop/dt)
dB <- rnorm(N,sd=sqrt(dt))

X0 <- 1

r <- 0.05
sigma <- 0.1

Tvec <- seq(0,Tstop,dt)
B <- c(0,cumsum(dB))

X <- X0*exp((r-0.5*sigma^2)*Tvec+sigma*B)

Tnow <- Tstop/2
Xnow <- X[N/2+1]

Tpred <- seq(Tstop/2,Tstop,dt)

Xpred <- Xnow*exp(r*(Tpred-Tnow))

Xlo <- Xnow*exp((r-0.5*sigma^2)*(Tpred-Tnow) - 1*sigma*sqrt(Tpred-Tnow))
Xhi <- Xnow*exp((r-0.5*sigma^2)*(Tpred-Tnow) + 1*sigma*sqrt(Tpred-Tnow))


plot(Tvec,X,type="l",xlab="Time",ylab="Price",ylim=c(0,max(c(Xhi,X))))
polygon(c(Tpred,rev(Tpred)),c(Xlo,rev(Xhi)),col="grey",border=NA)
lines(Tvec,X)
lines(Tpred,Xpred)



