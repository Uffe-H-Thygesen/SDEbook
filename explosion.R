pdf(file="explosion.pdf",width=6,height=4)

set.seed(12345)

Xmax <- 20
Xlim <- 10
N <- 3

T <- 10
dt <- 0.001

nt <- ceiling(T/dt)

tvec <- seq(0,T,length=nt+1)

dB <- array(rnorm(N*nt,sd=sqrt(dt)),c(nt,N))

B <- apply(dB,2,function(x)c(0,cumsum(x)))

X <- tan(B)

tau <- apply(X,2,function(x)sum(cumprod(abs(x)<Xmax)))

plot(c(0,max(tau))*dt,c(-1,1)*Xlim,type="n",xlab="t",ylab=expression(X[t]))

for(i in 1:N) lines(tvec[1:tau[i]],X[1:tau[i],i],lwd=2)

dev.off()
