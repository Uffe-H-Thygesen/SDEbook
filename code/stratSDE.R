### Sample BM

T <- 2*pi
N <- 64

### Sample random Fourier coefficients
V <- (rnorm(2*N+1) + 1i*rnorm(2*N+1))/(sqrt(2*pi))

### Associated frequencies
omega <- (-N):N

### Function to evaluate B.M. at a given time, for the V coefficients in globalenv
evalB <- function(t)
    {
        ## Hack to deal with 0 frequency
        omega <- omega + 1e-8

        return(sum(Re((exp(1i*omega*t)-1)*(V/(1i*omega)))))
    }

evalW <- function(t) sum(Re(exp(1i*omega*t)*V))

    
tvec <- seq(0,T,length=10001)
Bvec <- sapply(tvec,evalB)

dt <- diff(tvec[1:2])

plot(tvec,Bvec,type="l")

iBdB1 <- c(0,cumsum(Bvec[-1]*(diff(Bvec))))
iBdB2 <- c(0,cumsum(Bvec[-length(Bvec)]*(diff(Bvec))))
iBdB3 <- 0.5*Bvec^2

integrand <- function(tvec)sapply(tvec,function(t)evalB(t)*evalW(t))

iBdB4 <- integrate(integrand,lower=0,upper=T)

plot(tvec,iBdB1,type="l")
lines(tvec,iBdB2,col="red")
lines(tvec,iBdB3,col="blue")

