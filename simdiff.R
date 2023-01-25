## Simulating diffusion for "transport" chapter

set.seed(123456) # For reproducible results

Npart <- 100
xlim <- c(0,1)
ylim <- c(0,1)

## Diffusivity
D <- 1

## Simulation control
dt <- 0.00001
Nt <- 500

## Helper to project back into domain
project <- function(z,zlim)
  {
    I.left <- (z < zlim[1])
    z[I.left] <- 2*zlim[1]-z[I.left]

    I.right <- (z > zlim[2])
    z[I.right] <- 2*zlim[2]-z[I.right]

    return(z)
  }
    

## Init
xmean <- 0.5
ymean <- 0.5
xsd <- 0.05
ysd <- 0.05

X <- Y <- array(0,c(Npart,Nt))

X[,1] <- project(rnorm(Npart,xmean,xsd),xlim)
Y[,1] <- project(rnorm(Npart,ymean,ysd),ylim)

for(i in 2:Nt)
  {
    X[,i] <- project(X[,i-1] + sqrt(2*D*dt)*rnorm(Npart),xlim)
    Y[,i] <- project(Y[,i-1] + sqrt(2*D*dt)*rnorm(Npart),ylim)
  }

### Eulerian
Nf <- 50    # Number of Fourier coefficients

k <- 0:Nf

xvec <- seq(0,1,length=101)

FourierSol <- function(t)
  {
    phiF <- rep(1,length(xvec))
    for(i in k[-1])
      phiF <- phiF + 2*(-1)^i*exp(-D*(2*pi*i)^2*t)*cos(2*pi*i*xvec)

    return(phiF)
  }


### Eulerian densities

t0 <- xsd^2/2/D
phi0F <- FourierSol(t0)
phi0 <- outer(phi0F,phi0F)

maxPhi <- max(phi0F)^2

t1 <- t0+Nt*dt
phi1F <- FourierSol(t1)
phi1 <- outer(phi1F,phi1F)

### Create plot
require(rasterpdf) ## To avoid raster appearance (!)
raster_pdf(file="simdiff.pdf")

par(mfrow=c(2,2),mar=rep(0.1,4),cex=1.5)

## Color codes. White background for zero density:
col <- gray(seq(1,0,length=2^6+1))

trans <- function(x) x^0.4 # sqrt
zlim <- c(0,trans(maxPhi))

## Subpanels with Eulerian view
image(trans(abs(phi0)),zlim=zlim,col=col,xaxt="n",yaxt="n",bty="o") 
box()
text(xmean,ylim[2],"Initial",pos=1)
text(xlim[1],ymean,"Concentration",srt=90,adj=c(0.5,1.5))
image(trans(abs(phi1)),zlim=zlim,col=col,xaxt="n",yaxt="n",bty="o")
box()
text(xmean,ylim[2],"Terminal",pos=1)

## Subpanel with particles
cex <- 0.5 
plot(X[,1],Y[,1],pch=16,cex=cex,
     xlim=xlim, ylim=ylim,col="darkgrey",xaxt="n",yaxt="n")
text(xlim[1],ymean,"Particles",srt=90,adj=c(0.5,1.0))

# Pick a partcle which ends near desired end point 
m <- which.min((X[,Nt]-1)^2+(Y[,Nt]-1)^2)

plot(X[,Nt],Y[,Nt],pch=16,cex=cex,xlim=xlim,
     ylim=ylim,col="darkgrey",xaxt="n",yaxt="n")

points(X[m,],Y[m,],col="black",pch=16,cex=0.5*cex)

## points(X[m,1],Y[m,1],pch="+",cex=cex,col="green")
## points(X[m,Nt],Y[m,Nt],pch="+",cex=cex,col="red")

dev.off()
