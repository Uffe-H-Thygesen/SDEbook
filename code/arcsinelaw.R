### Simulate Brownian motion for t in [0,1] and return the time of last change of sign
simTau <- function(n=1,dt=1e-6,do.plot=FALSE)
{
    taus=numeric(n)
    for(i in 1:n)
        {
            dB=rnorm(1/dt)*sqrt(dt);
            B<- cumsum(dB);
            tau<- 1-dt*sum(cumprod(sign(rev(B))==sign(B[length(B)])));
            taus[i] <- tau
        }

    if(do.plot)
        {
            plot(seq(0,1,dt),c(0,B),type="l",main=tau)
            abline(h=0)
        }
    return(taus)
}

### Simulate a number of such last hitting times
taus <- simTau(n=1e3,dt=1e-3,do.plot=FALSE)

### Plot empirical c.d.f.
plot(sort(taus),seq(0,1,length=length(taus)),type="l")

### Compare with analytical result
plot(function(x)2/pi*asin(sqrt(x)),from=0,to=1,add=TRUE,col="red")

### Compare with (wrong?) analytical result
plot(function(x)2/pi*atan(sqrt(x)),from=0,to=1,add=TRUE,col="blue")
