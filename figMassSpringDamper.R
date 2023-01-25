### R-code to generate figures for the mass-spring-damper example
###
###

### Parameters
M <- 1          # Mass
K <- 1          # SPring constant
C <- 0.4        # Damper

P <- 1          # Impulse magnitude for impulse response
F0 <- 1         # Force apmplitude for frequency response

### Compute determinant
D <- C^2 - 4*M*K 

graphics.off()

if(D >= 0) stop("The system is overdamped ... I hadn't expected that to happen.")

lambda <- -C/2/M
omega <- sqrt (-D)/2/M
x <- function(t) P/M / omega * exp(lambda*t)*sin(omega*t)

Tmax <- 25/sqrt(lambda^2+omega^2)

pdf("impulse-response.pdf",width=3.5,height=3.5)
par(mar=c(2.6,3,1.5,0)+0.1,mgp=c(1.5,0.5,0))
plot(x,from=0,to=Tmax,xlab="Time [s]",ylab = "Position [m]",main="Impulse response")
dev.off()

omegas <- seq(omega/5,5*omega,omega/100)[-1]

pdf("frequency-response.pdf",width=3.5,height=3.5)
Hxs <- F0/(K+C*1i*omegas - M*omegas^2)
Hvs <- F0*1i*omegas/(K+C*1i*omegas - M*omegas^2)

par(mfrow = c(2,1),mar=c(1.1,3.1,1.6,0.2),mgp=c(1.5,0.5,0))
plot(omegas,Mod(Hxs),ylim=c(0.2,max(Mod(Hxs))),type="l",log="xy",xlab="",yaxt="n",xaxt="n",
     ylab="Amplitude [m]",main="Frequency response")
axis(side=2,at = c(0.25,1))

##  lines(omegas,Mod(Hvs),type="l",log="xy")
par(mar=c(2.6,3.1,0.1,0.2))
plot(omegas,Arg(Hxs),type="l",log="x",xlab="Angular frequency [rad/s]",
     ylab="Phase [rad]",yaxt="n",ylim=c(-pi,0))
axis(side=2,at = c(-pi,0),labels=c(expression(-pi),0))
##  lines(omegas,Arg(Hvs),type="l",log="xy")
dev.off()

###################################################
### Figure of the acf and spectrum of the force
pdf(file="force-acf-spectrum.pdf",width=8,height=4)

par(cex=3,mfrow=c(1,2))

### Range of plot
maxlag <- 5

### Difference at which time lags?
tauvec <- c(0.1,0.5)
sigmavec <- sqrt(1/tauvec)

acf <- function(tau,sigma,h) sigma^2*exp(-abs(h)/tau)

plot(function(h) acf(tauvec[1],sigmavec[1],h),from=-maxlag,to=maxlag,
     xlab="Time lag h",ylab=expression("A.c.f. "*rho[U](h)),lty=1,lwd=3)

for(i in 2:length(tauvec))
  plot(function(h) acf(tauvec[i],sigmavec[i],h),from=-maxlag,to=maxlag,add=TRUE,lty=i,lwd=3)


legend("topright",sapply(1:length(tauvec),function(i)as.expression(bquote(tau == .(tauvec[i])))),       lty=1:length(tauvec),lwd=3,horiz=FALSE)

## Spectrum
S <- function(w,sigma,tau) 2*sigma^2*tau/(1+w^2*tau^2)
maxfreq <- 10
minfreq <- 0.1

plot(function(w) S(w,sigmavec[1],tauvec[1]),from=minfreq,to=maxfreq,
     xlab=expression("Frequency "*omega),ylab=expression("Spectrum "*S[U](omega)),lty=1,lwd=3,log="xy",
     ylim=c(S(maxfreq,sigmavec[2],tauvec[2]),S(0,sigmavec[1],tauvec[1])))


for(i in 2:length(tauvec))
  plot(function(w) S(w,sigmavec[i],tauvec[i]),from=minfreq,to=maxfreq,add=TRUE,lty=i,lwd=3)

legend("bottomleft",sapply(tauvec,function(i)as.expression(bquote(tau == .(i)))),       lty=1:length(tauvec),lwd=3,horiz=FALSE)

dev.off()
###################################################


### Simulation; short run with long time between jumps 
set.seed(2345)
dt <- 0.01       # Time step
tau <- 15     # Expected time between jumps

rf <- rnorm     # Function for generating random force after a jump

Tstop <- 250
N <- ceiling(Tstop/dt)
Tstop <- N*dt
Tvec <- seq(0,Tstop,dt)

# Initialize arrays
X <- V <- F <- numeric(N+1)

F[1] <- rf(1)

for(i in 2:N+1)
  {
    ## Check if force has jumped
    if(rexp(1,rate=1/tau)<dt) F[i] <- rf(1) else F[i] <- F[i-1]

    V[i] <- V[i-1]*(1-C/M*dt) - K/M*X[i-1]*dt + F[i-1]*dt
    X[i] <- X[i-1]+(V[i]+V[i-1])/2*dt
  }

XVFts <- ts(data.frame(q=X,v=V,U=F),start=0,frequency=1/dt)


### Plot time series
pdf(file="figMassSpringDamper.pdf",height=3.5,width=6)
par(mar=c(4,4,0,0)+0.1)
plot(XVFts,main="")
dev.off()


### Simulation; long run with short time between jumps 
dt <- 0.01       # Time step
tau <- 0.1      # Expected time between jumps

rf <- rnorm     # Function for generating random force after a jump

Tstop <- 2500
N <- ceiling(Tstop/dt)
Tstop <- N*dt
Tvec <- seq(0,Tstop,dt)

# Initialize arrays
X <- V <- F <- numeric(N+1)

for(i in 2:N+1)
  {
    ## Check if force has jumped
    if(rexp(1,rate=1/tau)<dt) F[i] <- rf(1) else F[i] <- F[i-1]

    V[i] <- V[i-1]*(1-C/M*dt) - K/M*X[i-1]*dt + F[i-1]*dt
    X[i] <- X[i-1]+(V[i]+V[i-1])/2*dt
  }


XVFts <- ts(data.frame(q=X,v=V,F=F),start=0,frequency=1/dt)

XVFs <- spectrum(XVFts,plot=FALSE)
dat.array <- array(c(X,V,F),c(length(X),3))
spec <- mvfft(dat.array)
freqs <- seq(0,5,length=ceiling(dim(spec)[1]/2))

omegas <- seq(1/tau/100,2.5*max(c(omega,1/tau)),omega/100)[-1]
Hxs <- F0/(K+C*1i*omegas - M*omegas^2)
Hvs <- F0*1i*omegas/(K+C*1i*omegas - M*omegas^2)
Sff <- 2*tau/(1+omegas^2*tau^2)
Sxx <- Sff * Mod(Hxs)^2
Svv <- Sff * Mod(Hvs)^2

fact <- 2 # sqrt(2*pi) #dt /sqrt(2*pi)
# par(mfrow=c(3,1))

pdf(file="figMassSpringDamperSPEC.pdf",width=5,height=7)
layout(array(1:3,c(3,1)),width=1,heights=c(0.35,0.3,0.35))
par(mar=c(0,6,4,1),las=1)

plot(omegas,Sxx,type="l",log="xy",xaxt="n",
     xlab="",main="Variance spectra",
     ylab="")
title(ylab=expression("Position "*S[xx]*"[m"^2*s*"]"),line=4.5)
fs <- cut(log(XVFs$freq+0.1),100)
points(2*pi*exp(tapply(log(XVFs$freq),fs,mean)),fact*exp(tapply(log(XVFs$spec[,1]),fs,mean)),pch=16)
grid()

par(mar=c(0,6,1,1))
plot(omegas,Svv,type="l",log="xy",xlab="",xaxt="n",
     ylab="")
title(ylab=expression("Velocity "*S[vv]*"[m"^2/s*"]"),line=4.5)
points(2*pi*exp(tapply(log(XVFs$freq),fs,mean)),fact*exp(tapply(log(XVFs$spec[,2]),fs,mean)),pch=16)
grid()

par(mar=c(5,6,1,1))
if(tau<1) ylims <- Sff[1]*c(0.1,10) else ylims <- NULL
plot(omegas,Sff,type="l",log="xy",
     xlab="Frequency [rad/s]",ylab="",ylim=ylims)
title(ylab=expression("Force "*S[UU]*"[m"^2/s^3*"]"),line=4.5)
points(2*pi*exp(tapply(log(XVFs$freq),fs,mean)),fact*exp(tapply(log(XVFs$spec[,3]),fs,mean)),pch=16)
grid()
dev.off()
