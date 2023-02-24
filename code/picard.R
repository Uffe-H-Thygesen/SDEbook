lambda <- 1
sigma <- 0.5

f <- function(x) lambda * x
g <- function(x) sigma * x

StratInt <- function(F,M)
    c(0,cumsum((F[-1]-0.5*diff(F))*diff(M)))

x0 <- 1

T <- 3
dt <- 1e-3
nt <- round(T/dt)

t <- seq(0,T,length=nt+1)

dB <- rnorm(nt,sd=sqrt(dt))
B <- c(0,cumsum(dB))

# Cheat and construct a Brownian bridge, in order to make sure the plot comes out nice
B <- B - B[nt+1]*t/T
dB <- diff(B)

Xa <- x0*exp(lambda*t+sigma*B)

N <- 4

Xi <- array(0,c(nt+1,N+1))

Xi[,1] <- x0
for(i in 1:N)
    Xi[,i+1] <- x0 + StratInt(f(Xi[,i]),t) + StratInt(g(Xi[,i]),B)

pdf(file="picard.pdf",width=8)

plot(1.1*range(t),range(c(range(Xi),range(Xa))),type="n",xlab="t",ylab="x")

for(i in 1:(N+1))
{
    lines(t,Xi[,i],lty=1,col="grey")
    text(T,Xi[nt,i],pos=4,substitute(X[t]^k,list(k=i-1)))
}

lines(t,Xa,lwd=1)

text(T,Xa[nt],pos=4,expression(X[t]))

dev.off()
