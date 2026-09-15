## The stochastic van der Pol oscillator

require(SDEtools)
set.seed(123456)

## Define the model, each element in the vector equation
fx <- function(x,v) v
fv <- function(x,v) v*(mu - x^2)-x
gv <- function(x,v) sigma

## Construct the vector fields
f <- function(xv) c(fx(xv[1],xv[2]),fv(xv[1],xv[2]))
g <- function(xv) c(0,gv(xv[1],xv[2]))

## Parameters
sigma <- 0.5
mu <- 1

## Simulation control
T <- 100
dt <- 0.01
times <- seq(0,T,dt)

B <- rBM(times)

xv0 <- c(0,0)

sim <- euler(f,g,times,xv0,B)

## Plot trajectory and phase plane
pdf(file="vanderPol.pdf",width=6,height=4)
par(mfrow=c(1,2))
plot(times,sim$X[,1],type="l",xlab="Time t",ylab="Position X")
plot(sim$X[,1],sim$X[,2],type="l",xlab="Position X",ylab="Velocity V")
dev.off()

## Computing the stationary density

xmax <- 1.5*max(abs(sim$X[,1]))
vmax <- 1.5*max(abs(sim$X[,2]))

xgrid <- seq(-xmax,xmax,length=101)
vgrid <- seq(-vmax,vmax,length=91)

xc <- cell.centers(xgrid)
vc <- cell.centers(vgrid)

## Diffusivity in each direction
Dx <- function(x,v) 0
Dv <- function(x,v) 0.5*gv(x,v)^2

## Construct discretized generator
G <- fvade2d(fx,fv,Dx,Dv,xgrid,vgrid)

## Find the stationary distribution
pi <- unpack.field(StationaryDistribution(G),length(xc),length(vc))

## Tweaked image for raster-free visualization
my.image <- function(x,y,z,col=c(0,1),xlab="",ylab="",main="")
{
    plot(range(x),range(y),type="n",xlab=xlab,ylab=ylab,
         main=main,xaxs="i",yaxs="i")
    rangex <- range(x)
    rangey <- range(y)
    rangez <- range(z)
    scale <- function(z) (z-rangez[1])/diff(rangez)*diff(col)+col[1]
    zraster <- 
    rasterImage(t(scale(z[nrow(z):1,])),rangex[1],rangey[1],rangex[2],rangey[2],col=col)
}

pdf(file="vanderPol-density.pdf",width=5,height=4)
my.image(xgrid,vgrid,t(pmax(pi,0))^0.75,col=c(1,0),
         xlab="Position X",ylab="Velocity Y")
dev.off()

