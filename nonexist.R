# Figure to demonstrate explosions

par(mfrow=c(1,2),cex=1.25)

# DETERMINISTIC

plot(function(t)tan(t),ylim=c(-10,10),from=0,to=3/2,lwd=2,xlab="Time t",ylab="x")


# Stochastic
N <- 10

T <- 5
dt <- 0.01
Nt <- ceiling(T/dt-1e-5)

t <- seq(0,T,length=Nt+1)

dB <- array(sqrt(dt)*rnorm(Nt*N),c(Nt,N))
B <- apply(dB,2,cumsum)

B <- rbind(rep(0,N),B)

In  <- (abs(B)<pi/2)
In <- apply(In,2,cumprod)

X <- tan(B*In)
X[!In] <- NA

UseThese <- c(1,6,9)

AnyIn <- apply(In[,UseThese],1,function(x)as.numeric(any(as.logical(x))))


Tmax <- sum(AnyIn)*dt

plot(c(0,Tmax),c(-10,10),type="n",xlab="Time t",ylab="x")
for(i in 1:length(UseThese))
    lines(t,X[,UseThese[i]],col=i)

