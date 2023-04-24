## Simulation of the noisy harmonic oscillator
require(SDEtools)
require(Matrix)

set.seed(1234567)

omega <- 1  ## Natural frequency
x0 <- 1     ## Initial condition

## Simulation time control
T <- 100
dt <- 0.1
tv <- seq(0,T,dt)

B <- rvBM(tv,n=2)

pdf(file="harmonic.pdf",width=5,height=7)

par(mfcol=c(3,2),mar=c(5,5,0.1,0.1),mgp=c(2,1,0))

## Repeat for different values of the damping
for(lambda in c(0.05,0.5))
{
    ## Scale noise intensity to get the same stationary variance
    sigma <- sqrt(2*lambda)

    ## System matrix
    A <- array(c(-lambda,omega,-omega,-lambda),c(2,2))
    I <- diag(c(1,1))

    f <- function(x)A %*% x
    g <- function(x)sigma*I

    sim <- euler(f,g,tv,c(0,0),B=B)

    plot(tv,sim$X[,1],type="l",
         xlab="Time t",ylab=expression(X[t]^(1)),
         ylim=4*c(-1,1))


    ## Tabulate and plot the a.c.f.
    acf <- function(h)
    {
        if(h >= 0)
            return( sigma^2/2/lambda * as.array(expm(A*h)) ) else 
            return(t(acf(-h)))
    }

    hmax <- 20
    hv <- seq(-hmax,hmax,dt)
    acftab <- array(NA,c(length(hv),2,2))
    for(i in 1:length(hv)) acftab[i,,] <- acf(hv[i])

    plot(hv,acftab[,1,1],type="l",ylim=c(-1,1),
         xlab="Time lag h",ylab=expression("A.c.f. "*rho[11](h)))

    ## Tabulate and plot the spectrum
    spec <- function(w)
        return(sigma^2* solve(-1i*w*I - A) %*% t(solve(1i*w*I-A)))

    wv <- omega*exp(seq(-2,2,length=1001))
    
    st1 <- spectab <- array(NA,c(length(wv),2,2))

    for(i in 1:length(wv)) spectab[i,,] <- spec(wv[i])

    plot(wv,Re(spectab[,1,1]),type="l",log="xy",
         xlab=expression("Frequency "*omega),ylab=expression("Spectrum "*S[11](omega)),
         ylim=c(0.01,20))


}

dev.off()
