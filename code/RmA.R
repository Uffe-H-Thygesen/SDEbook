require(SDEtools)

graphics.off()

Nmax <- 1.2
Pmax <- 0.8

nn <- 151
np <- 152

ni <- seq(0,Nmax,length=nn)
pi <- seq(0,Pmax,length=np)

nc <- 0.5*(head(ni,-1)+tail(ni,-1))
pc <- 0.5*(head(pi,-1)+tail(pi,-1))

## Data: Drift from a Rosenzwerig-MacArthur model.
## Diffusivity corresponding to multiplicative noise
r <- 1.
K <- 1.
Cmax <- 1

epsilon <- 0.3
mu <- 0.15 
beta <- 4

Dn <- function(n,p) 0.001*n^2
Dp <- function(n,p) 0.001*p^2

fn <- function(n,p) r*n*(1-n/K)-Cmax*n*p/(n+Cmax/beta)
fp <- function(n,p) epsilon*Cmax*n*p/(n+Cmax/beta) - mu*p

G <- fvade2d(ux=function(n,p) fn(n,p) - 2*Dn(1,1)*n,
             uy=function(n,p) fp(n,p) - 2*Dp(1,1)*p,
             Dx=Dn,
             Dy=Dp,
             ni,pi)

rho <- StationaryDistribution(G)

rho <- unpack.field(rho,length(nc),length(pc))

## Sim
times <- seq(0,1e3,0.1)

sn <- sqrt(Dn(1,1)*2)
sp <- sqrt(Dp(1,1)*2)
    
sim <- euler(f=function(np) c(fn(np[1],np[2]),fp(np[1],np[2])),
             g=function(np) diag(np*c(sn,sp)),
             times=times,
             x0=c(0.01,0.01),
             p=abs)

nt <- 1000

my.image <- function(x,y,z,col=c(1,0),xlab="",ylab="",main="")
{
    plot(ylim,ylim,type="n",xlab=xlab,ylab=ylab,
         main=main,asp=1,ylim=ylim,xlim=ylim)
    rangex <- range(x)
    rangey <- range(y)
    rangez <- range(z)
    scale <- function(z) (z-rangez[1])/diff(rangez)*diff(col)+col[1]
    zraster <- 
    rasterImage((scale(z[nrow(z):1,])),rangex[1],rangey[1],rangex[2],rangey[2],col=col)
}

ylim <- c(0,max(nc))

pdf(file=paste("epsilon-",epsilon,".pdf",sep=""),height=2,width=6)
par(mfrow=c(1,3),mar=c(2,2,0.5,0.5),tck=0.025,mgp=c(1,0.1,0))

my.image(nc,pc,rho,
      xlab="Prey",ylab="Predators")

plot(tail(sim$X[,1],3*nt),tail(sim$X[,2],3*nt),type="l",asp=1,
     xlim=ylim,ylim=ylim,xlab="Prey",ylab="Predators")
matplot(tail(times,nt),tail(sim$X,nt),type="l",lty=1:2,col=1,ylim=ylim,
        xlab="Time",ylab="Prey, Predators")
legend("topright",lty=1:2,legend=c("Prey","Predators"))
dev.off()

## dev.new()
## par(mfrow=c(1,2))
## hist(sim$X[,1],freq=FALSE)
## lines(nc,apply(rho,2,sum)/diff(ni))

## hist(sim$X[,2],freq=FALSE)
## lines(pc,apply(rho,1,sum)/diff(pi))
