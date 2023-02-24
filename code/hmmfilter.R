## Figure for the notes about state estimation in the Cox-Ingersoll-Ross process.
## Adaptation of exercise 8 from course 02425

## Start from a clean slate
rm(list=ls())
graphics.off()

set.seed(12345)

require(Matrix)
require(fields)
require(SDEtools)

lambda <- 1
xi <- 1
gamma <- 1

# Define model. Note abs to handle negative x.
f = function(x) lambda*(xi-x);
g = function(x) gamma*sqrt(abs(x));

# Advection-diffusion form
D = function(x) 0.5*gamma^2*x;
u = function(x) f(x) - 0.5*gamma^2;

#############################################################################
## Grid
#############################################################################

## Use a grid which has finer resolution near x=0 (because we can)
xii = seq(0,2,0.025)^2                                # Cell Interfaces
dx <- diff(xii)
xc = 0.5*(tail(xii,-1) + head(xii,-1))   # Cell centers


#############################################################################
## Simulate using Euler method
#############################################################################

dt = 0.01;
Tmax = 20;

tvec = seq(0,Tmax,dt)
nt <- length(tvec)

# Euler scheme for simulation, with reflection at 0 to avoid negative x's
B <- rBM(tvec)
sim <- euler(f,g,tvec,xi,B,abs)

x <- sim$X

## Export states
sim <- data.frame(t=tvec,X=x)
write.table(sim,"hmm-states.txt",quote=FALSE,row.names=FALSE)

## Generate random measurements
tsample = 0.1;
vsample = 1;

sampleIndeces = round(seq(1,length(tvec),(tsample/dt)))
tm = tvec[sampleIndeces];
ymean = vsample*x[sampleIndeces];

# Generate random measurements 
U <- runif(length(tm))
ymeas <- qpois(p=U,lambda=ymean)

## Export measurements
obs <- data.frame(t=tm,Y=ymeas)
write.table(obs,file="hmm-obs.txt",quote=FALSE,row.names=FALSE)

#############################################################################
## HMM filter
#############################################################################

# Data likelihood
dl = function(x,y) dpois(y,lambda=vsample*x)

fit <- HMMfilterSDE(u,D,xii,'r',NULL,tm,ymeas,function(x,y)dpois(y,lambda=x*vsample),do.smooth=TRUE)

## Show confidence intervals, mean state, and true state 

ci <- 1/6
psiCDF <- t(apply(fit$psi,1,cumsum))

Xl <- apply(psiCDF,1,function(cdf)approx(cdf,xc,ci)$y)
Xu <- apply(psiCDF,1,function(cdf)approx(cdf,xc,1-ci)$y)

## Compute and posterior mean
XestMean <- fit$psi %*% xc

## Illustrate one time step: Find one with small mean and a positive measurement
o <- order(XestMean)
ii <- o[(which(ymeas[o] > 0))[2]]

Xsmooth <- apply(fit$pi,1,function(p)sum(p*xc))    

piCDF <- t(apply(fit$pi,1,cumsum))

Xsl <- xc[1+apply(piCDF,1,function(x)sum(x<ci))]
Xsu <- xc[1+apply(piCDF,1,function(x)sum(x< (1-ci)))]


#############################################################################
## Likelihood inference
#############################################################################

## Define likelihood function of parameter

loglik <- function(xi)
{
    ## Define model. Note abs to handle negative x.
    f = function(x) lambda*(xi-x);
    g = function(x) gamma*sqrt(abs(x));
    
    ## Advection-diffusion form
    D = function(x) 0.5*gamma^2*x;
    u = function(x) f(x) - 0.5*gamma^2;

    fit <- HMMfilterSDE(u,D,xii,'r',
                        function(x)x,tm,ymeas,function(x,y)dpois(y,lambda=x*vsample),
                        do.smooth=TRUE,
                        pfun=function(x,y)ppois(y,vsample*x))
    return(sum(log(fit$c)))
}

xis <- seq(0.25,2,0.05)
logliks <- sapply(xis,loglik)
est.xi <- xis[which.max(logliks)]

est.xi <- nlminb(est.xi,function(p)-loglik(p),lower=0,upper=2)$par

Il <- xis < est.xi
Iu <- xis > est.xi

CLu <- approx(logliks[Iu],xis[Iu],max(logliks) - 0.5*qchisq(0.95,df=1))$y
CLl <- approx(logliks[Il],xis[Il],max(logliks) - 0.5*qchisq(0.95,df=1))$y

conf.int <- c(CLl,CLu)


#############################################################################
## TMB Model
#############################################################################

require(TMB)
compile("cir_tmb.cpp")

dyn.load(dynlib("cir_tmb"))

# Data for TMB
data <- list(tsim=tvec,iobs=sampleIndeces-1,Y=ymeas)

# Initial geuess on latent variables
Z <- sqrt(x)

parameters <- list(
    Z=Z,
    lambda=lambda,
    xi=xi,
    loggamma=log(gamma),
    logv=log(vsample)
    )

fixed <- factor(NA)

obj <- MakeADFun(data,parameters,map=list(logv=fixed),random=c("Z"),DLL="cir_tmb")

# Estimate latent variables for true parameters
if(TRUE)
{
    obj$fn()
    sdr <- sdreport(obj)

    ## Generate one step predictions
    pred  <- oneStepPredict(obj,observation.name="Y",data.term.indicator="keep",
                            discrete=TRUE,method="oneStepGeneric",range=c(0,Inf),
                            parallel=TRUE)
}

# Get predictions of states with std.dev.
predTP <- summary(sdr,"random")
ZpredTP <- predTP[,1]
ZsdTP <- predTP[,2]

system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))
sdr <- sdreport(obj)

################################################################################################################
## Plots


## Data
pdf(file="hmm-data.pdf",width=6,height=4)
par(mar = c(5, 4, 4, 4) + 0.1)
xmax <- 5
plot(tvec,x,type="l",xlab="Time",ylab='Abundance',ylim=c(0,xmax))
par(new=TRUE)
plot(tm,ymeas,pch='o',ylim=c(0,xmax*vsample),axes=FALSE,ann=FALSE)
axis(4)
box()
mtext("Observations", side=4, line=3)
dev.off()

## Schematic of the algorithm
pdf(file="hmm-updates.pdf",width=8,height=4)
par(mfcol=c(1,2),mar=c(5,4,4,4)+0.1)
plot(xc,fit$psi[ii-1,]/dx,type="l",xlab="x",ylab="P.d.f.",ylim=c(0,1.7))
lines(xc,fit$phi[ii,]/dx,lty="dashed")
legend("topright",legend=c(expression(psi[i]),expression(phi[i+1])),lty=c("solid","dashed"))

plot(xc,fit$phi[ii,]/dx,type="l",xlab="x",ylab="P.d.f.",ylim=c(0,1.7))
lines(xc,fit$psi[ii,]/dx,lty="dashed")
par(new=TRUE)
plot(xc,dl(xc,ymeas[ii]),axes=FALSE,ann=FALSE,type="l",lty="dotted")
axis(4)
box()
mtext("State likelihood", side=4, line=3)

legend("topright",legend=c(expression(phi[i+1]),expression(l[i+1]),expression(psi[i+1])),lty=c("solid","dotted","dashed"))
dev.off()


## HMM estimate 
pdf(file="hmmfilter.pdf",width=6,height=4)
plot(range(tvec),range(c(0,Xu,x)),type="n",xlab="t",ylab="x")
polygon(c(tm,rev(tm)),c(Xl,rev(Xu)),col="grey",border=NA)
lines(tm,XestMean)
lines(tvec,x,type="l",lwd=1,lty="dashed")
legend("topright",legend=c("True state","Estimated state"),lty=c("dashed","solid"),lwd=c(1,1))
dev.off()

## Smoothing estimate obtained with HMM
pdf(file="hmm-smooth.pdf",width=6,height=4)
plot(range(tvec),range(c(0,Xu,x),na.rm=TRUE),type="n",xlab="t",ylab="X")
polygon(c(tm,rev(tm)),c(Xsl,rev(Xsu)),col="grey",border=NA)
lines(tm,XestMean,lty="solid")
lines(tm,Xsmooth,lty="dashed",lwd=2)
legend("topright",legend=c("Smoothed state","Estimated state"),lty=c("dashed","solid"),lwd=c(2,1))
dev.off()

## Likelihood estimation
pdf(file="hmm-loglik.pdf",width=6,height=4)
plot(xis,logliks,type="l",xlab=expression(xi),ylab=expression(log*Lambda(xi)))
polygon(rep(conf.int,c(2,2)),(range(logliks)+c(-5,5))[c(1,2,2,1)],col="grey",border=NA)
lines(xis,logliks)
abline(v=xi,lty="dashed")
points(est.xi,max(logliks),pc=16,cex=3)
dev.off()

## Compare HMM with TMB
pdf(file="hmm-tmb.pdf",width=6,height=4)
plot(tvec,ZpredTP^2,type="l",lwd=2,xlab="t",ylab="X")
lines(tm,Xsmooth)

## lines(tm,apply(fit$pi,1,function(p)xc[which.max(p)]),lty="dashed")
## lines(tm,xc[1+apply(piCDF,1,function(x)sum(x<0.5))],lty="dotted")

legend("top",legend=c("TMB","HMM, mean"),lty=c("solid","solid"),lwd=c(2,1),ncol=1)

dev.off()

