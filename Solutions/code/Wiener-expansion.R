### Demonstration of Wiener expansion of Brownian motion
###
### i.e., sampling BM in terms of Fourier coefficients,
### first simulating white noise, then integrating.

### Max no. frequencies
N <- 64

### Number of time steps for plotting
nt <- 1e2
dt <- 2*pi/nt

tvec <- seq(0,2*pi,length=nt+1)
tvec <- c(0.5,1,1.5,2,2.5,5)

omega <- (-N):N

simB <- function()
    {
        V <- (rnorm(2*N+1) + 1i*rnorm(2*N+1))/(sqrt(2*pi))
        expIwt <- exp(outer(1i*omega,tvec))
        Wn <- expIwt*V

        W <- V %*% expIwt
        ## plot(tvec,Re(W),type="l")
        
        ## Analytical integration with a hack
        omega <- omega + 1e-8
        IexpIwt <- (exp(outer(1i*omega,tvec))-1)/(1i*omega)
        B <- V %*% IexpIwt
        B <- as.numeric(Re(B)) + 1i * as.numeric(Im(B))
        return(B)
    }

M <- 1000
Bs <- sapply(1:M,function(i)simB())

graphics.off()

plot(tvec,Re(Bs[,1]),type="l",ylim=3*sqrt(2*pi)*c(-1,1))
for(i in 1:M) lines(tvec,Re(Bs[,i]),col="red")
for(i in 1:M) lines(tvec,Re(Bs[,i]),col="black")

dev.new()
EReB <- apply(Re(Bs),1,mean)
EImB <- apply(Im(Bs),1,mean)
VReB <- apply(Re(Bs),1,var)
VImB <- apply(Im(Bs),1,var)

par(mfrow=c(2,1))
plot(tvec,EReB)
lines(tvec,EImB)
plot(tvec,VReB)
lines(tvec,VImB)

