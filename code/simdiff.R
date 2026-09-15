## Simulating diffusion for "transport" chapter

require(SDEtools)
require(Matrix)

set.seed(123456) # For reproducible results

Npart <- 100     ## Number of particles
xylim <- c(0,1)  ## Plotting region (x/y-symmetric)

## Diffusivity
D <- 2

## Simulation control
dt <- 1e-5
T <- 5e-3
tv <- seq(0,T,dt)
Nt <- length(tv)

## Init
xymean <- 0.5
xysd <- 0.05

XY0 <- rnorm(2*Npart,xymean,xysd)

XY <- rvBM(tv,n=2*Npart,sigma=sqrt(2*D),B0=XY0)
X <- XY[,1:Npart]
Y <- XY[,-(1:Npart)]

### Eulerian
xvec <- seq(0,1,length=101)

phi0 <- dnorm(xvec,xymean,xysd)  ## P.d.f. of each coordinate
PHI0 <- outer(phi0,phi0)         ## Joint p.d.f. of both coordinates

phiT <- dnorm(xvec,xymean,sqrt(xysd^2+2*D*T))
PHIT <- outer(phiT,phiT)         ## ... same for the terminal position

maxPhi <- max(phi0)^2            ## For scaling of the colors

### Create plot
require(rasterpdf) ## To avoid raster appearance (!)
raster_pdf(file="simdiff.pdf")

par(mfrow=c(2,2),mar=rep(0.1,4),cex=1.5)

## Color codes. White background for zero density. Increase contrast at low
## densities
col <- gray(seq(1,0,length=2^8+1)^4)

zlim <- c(0,maxPhi)

## Subpanels with Eulerian view
image(PHI0,zlim=zlim,col=col,xaxt="n",yaxt="n",bty="o") 
box()
text(xymean,xylim[2],"Initial",pos=1)
text(xylim[1],xymean,"Concentration",srt=90,adj=c(0.5,1.5))
image(PHIT,zlim=zlim,col=col,xaxt="n",yaxt="n",bty="o")
box()
text(xymean,xylim[2],"Terminal",pos=1)

## Subpanel with particles
cex <- 0.5 
plot(X[1,],Y[1,],pch=16,cex=cex,
     xlim=xylim, ylim=ylim,col="darkgrey",xaxt="n",yaxt="n")
text(xylim[1],xymean,"Particles",srt=90,adj=c(0.5,1.0))

# Pick a partcle which ends near desired end point 
m <- which.min((X[Nt,]-1)^2+(Y[Nt,]-1)^2)

plot(X[Nt,],Y[Nt,],pch=16,cex=cex,xlim=xylim,
     ylim=ylim,col="darkgrey",xaxt="n",yaxt="n")

lines(X[,m],Y[,m],col="black",pch=16,cex=0.5*cex)

dev.off()
