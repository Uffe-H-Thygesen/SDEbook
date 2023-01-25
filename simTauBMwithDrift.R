### Simulate the hitting time for BM with drift, X = u*t + B

u <- 1
x <- 1
N <- 10000
T <- 2

rtau <- function(u,x,dt=1e-3,T=2)
    {
        t <- seq(0,T,dt)
        dB <- rnorm(length(t)-1,0,sqrt(dt))
        B <- c(0,cumsum(dB))
        X <- u*t + B
        tau <- (sum(cumprod(X<x))-1)*dt
        return(tau)
    }


tau <- sapply(1:N,function(i)rtau(u=u,x=x,T=T))

tau[tau==T] <- Inf

hist(tau,freq=FALSE,xlab=expression(t),ylab="Density",col="lightgrey",
    main=NULL) # expression("Histogram of "*tau*" and theoretical p.d.f. "*f[tau](t)))

plot(function(t)x*t^(-3/2)*dnorm(x/sqrt(t) ) *exp( u* x - u^2*t/2),
     from=0,to=2,add=TRUE,lwd=3)

dev.copy2pdf(file="simTauBMwithDrift.pdf")
