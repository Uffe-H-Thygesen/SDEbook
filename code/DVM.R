## Diel vertical migration as an optimal control problem
graphics.off()
rm(list=ls())

require(SDEtools)
require(fields)

## Growth
h <- function(y) hmax/(1+exp((y-yNC)/wNC))

## Mortality at the surface (following light intensity)
I <- function(t) mu1 /(1 + exp(A*cos(t/T*2*pi)))

## Mortality
mu <- function(y,t) I(t) * exp(-kl*y) + mu0

## Parameters
hmax <- 1   # J/day ; energy harvest rate at surface
wNC <- 10   # m, width of nutricline
yNC <- 50   # m, position of nutricline
kl <- 0.05  # 1/m , absorption coefficient of light
A <- 13     # 1, day/night difference between light
T <- 1      # day; period
mu1 <- 0.5     # 1/day ; max predation mortality
mu0 <- 1e-2 # 1/day ; base predation mortality
Dy <- 500   # m^2/day ; diffusivity
nu <- 2e-5     # J*day/m^2 ; Cost of motion 

## Discretization ; number of cells
nt <- 120
ny <- 100

tv <- seq(0,1,length=nt+1)
tc <- cell.centers(tv)
dt <- diff(tv)

yv <- seq(0,2*yNC,length=ny+1)
yc <- cell.centers(yv)
dy <- diff(yv)

## Generator without control: Pure diffusion in vertical and time flies with speed 1
G0 <- fvade2d(ux=function(x,y)1,uy=function(x,y)0, Dx=function(x,y)0,Dy=function(x,y)Dy,
              xgrid=tv,ygrid=yv,
              bc=list(N="r",S="r",E="p",W="p"))
## Generators corresponding to control: Swimming down 
G1 <- fvade2d(ux=function(x,y)0,uy=function(x,y)1, Dx=function(x,y)0,Dy=function(x,y)0,
              xgrid=tv,ygrid=yv,
              bc=list(N="r",S="r",E="p",W="p"))
## Generators corresponding to control: Swimming up 
G2 <- fvade2d(ux=function(x,y)0,uy=function(x,y)-1,Dx=function(x,y)0,Dy=function(x,y)0,
              xgrid=tv,ygrid=yv,
              bc=list(N="r",S="r",E="p",W="p"))

## Mortalities 
MU <- outer(yc,tc,mu) 
M <- Diagonal(nt*ny,as.numeric(MU)) 

G0 <- G0 - M

## Running reward: Harvest minus cost of motion
h0 <- rep(h(yc),nt)
k <- function(u) h0 - 0.5*nu*rowSums(u^2)

## Optimal swimming speed
uopt <- function(dV) pmax(dV/nu,0)

sol <- PolicyIterationRegular(G0,G1=list(G1,G2),k,list(uopt,uopt),
                              do.minimize=FALSE,
                              do.return.QSD=TRUE)

## Fitness for an arbitrary strategy
Dy <- 5000

## Generator without control: Pure diffusion in vertical and time flies with speed 1
G0 <- fvade2d(ux=function(x,y)1,uy=function(x,y)0, Dx=function(x,y)0,Dy=function(x,y)Dy,
              xgrid=tv,ygrid=yv,
              bc=list(N="r",S="r",E="p",W="p")) - M 
u <- function(y,t) 0
U <- outer(yc,tc,Vectorize(u))
Gcl <- G0 + Diagonal(nrow(G0),pmax(0,U))*G1 + Diagonal(nrow(G0),pmax(0,-U))*G2
V0 <- -as.numeric(solve(Gcl,h0-0.5*nu*as.numeric(U^2)))

qsd.BM <- QuasiStationaryDistribution(Gcl)


my.image <- function(x,y,z,col=col,xlab="",ylab="",main="")
{
    plot(range(x),range(y),type="n",xlab=xlab,ylab=ylab,
         main=main,xaxs="i",yaxs="i",ylim=rev(range(y)))
    rangex <- range(x)
    rangey <- range(y)
    rangez <- range(z)
    scale <- function(z) (z-rangez[1])/diff(rangez)*diff(col)+col[1]
    zraster <- 
    rasterImage(t(scale(z[,ncol(z):1])),rangex[1],rangey[1],rangex[2],rangey[2],col=col)
}

my.plot <- function(Z,col = c(0.3,1),xlab="",ylab="",main="")
{
    if(is.null(dim(Z))) Z <- unpack.field(Z,nt,ny)
    my.image(tc,yc,t(Z),col=col,xlab=xlab,ylab=ylab,main=main)
    contour(tc,yc,t(Z),nlevels=5,add=TRUE)
}

## Graphical parameters
mgp <- c(1.5,0.45,0)
tcl <- -0.4
mar <- c(2.5,3,2.5,0.5)

pdf(file="DVM-Model.pdf",width=6,height=3)
par(mfrow=c(1,2),mgp=mgp, tcl=tcl, mar=mar)
plot(h(yc),yc,ylim=rev(range(yv)),type="l",xlab="Harvest h [J/day]",ylab="Depth x [m]",xaxs="i",yaxs="i")
my.plot(MU,main="Mortality [1/day]",xlab="Time t [day]",ylab="Depth x [m]")
dev.off()

pdf(file="DVM-BM.pdf",width=6,height=5)
par(mfrow=c(2,1),mgp=mgp, tcl=tcl, mar=mar)
my.plot(V0,main="Fitness V [J]",xlab="Time t [day]",ylab="Depth [m]")
my.plot(pmax(0,qsd.BM$vector)/dy*nt,main="Quasi-stationary density [1/m] ",xlab="Time [day]",ylab="Depth [m]",col=c(1,0.3))
dev.off()

pdf(file="DVM.pdf",width=6,height=9)
par(mfrow=c(3,1),mgp=mgp, tcl=tcl, mar=mar)
VV <- unpack.field(sol$V,nt,ny)
UU <- unpack.field(sol$u %*% c(1,-1),nt,ny)

my.plot(VV,main="Fitness V [J]",xlab="Time t [day]",ylab="Depth x [m]")
my.plot(UU,main="Downward speed [m/day]",xlab="Time [day]",ylab="Depth [m]")
my.plot(pmax(0,sol$qsd.vector)/dy*nt,main="Quasi-stationary density [1/m] ",xlab="Time [day]",ylab="Depth [m]",col=c(1,0.3))

dev.off()
