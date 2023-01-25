### Figure of the acf and spectrum of he difference quotient of BM

graphics.off()

pdf(file="BM-diff-quotient.pdf",width=8,height=4)

par(cex=3,mfrow=c(1,2))

### Range of plot
maxlag <- 2

### Difference at which time lags?
kvec <- c(0.1,1)

acf <- function(h,k) pmax((k-abs(h))/k^2 ,0)

plot(function(h) acf(h,kvec[1]),from=-maxlag,to=maxlag,
     ylim=c(0,acf(0,kvec[1])*1.1),
     xlab="Time lag h",ylab=expression("A.c.f. "*rho[X](h)),lty=1,lwd=3)

for(i in 2:length(kvec))
  plot(function(h) acf(h,kvec[i]),from=-maxlag,to=maxlag,add=TRUE,lty=i,lwd=3)

legend("right",legend=paste("k=",kvec),lty=1:length(kvec),lwd=3,horiz=FALSE)

## Spectrum
S <- function(w,k) 2*(1-cos(w*k))/k^2/w^2
maxfreq <- 19

plot(function(w) S(w,kvec[1]),from=-maxfreq,to=maxfreq,
     xlab=expression("Frequency "*omega),ylab=expression("Spectrum "*S[X](omega)),lty=1,lwd=3,
     ylim = c(0,1.1*S(1e-5,1)))

for(i in 2:length(kvec))
  plot(function(w) S(w,kvec[i]),from=-maxfreq,to=maxfreq,add=TRUE,lty=i,lwd=3)

legend("right",legend=paste("k=",kvec),lty=1:length(kvec),lwd=3,horiz=FALSE)

dev.off()
    
