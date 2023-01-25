### Figure "potentials" for chapter "potential"

graphics.off()

par(mfrow=c(2,2),mar=0.5*c(4, 4, 4, 2) + 0.1)

dx = 0.01
xv <- seq(-2,2,dx)

xr <- range(xv)
ur <- c(-0.5,2)

# Constant potential
plot(xr,ur,type="n",xlab="x",ylab="U(x)",main=expression(U[0](x)))
lines(c(-1,1),c(0,0),lwd=3)
lines(c(-1,-1),ur,lty="dashed")
lines(c(1,1),ur,lty="dashed")
text(-1.5,mean(ur),expression(U==Inf))
text(+1.5,mean(ur),expression(U==Inf))

text(min(xr)+0.1*diff(xr),max(ur),"a",cex=1.5)


# Linear potential
plot(xr,ur,type="n",xlab="x",ylab="U(x)",main=expression(U[1](x)))
lines(c(0,max(xr)),c(0,max(xr)),lwd=3)
lines(c(0,0),ur,lty="dashed")

text(-1.5,mean(ur),expression(U==Inf))

text(min(xr)+0.1*diff(xr),max(ur),"b",cex=1.5)

## Quadratic potential
plot(xr,ur,type="n",xlab="x",ylab="U(x)",main=expression(U[2](x)))
lines(xv,0.5*xv^2,lwd=3)

text(min(xr)+0.1*diff(xr),max(ur),"c",cex=1.5)


# Double-well
plot(xr,ur,type="n",xlab="x",ylab="U(x)",main=expression(U[4](x)))
U4 <- 0.25*xv^4 - 0.5*xv^2
lines(xv,U4,lwd=3)

text(min(xr)+0.1*diff(xr),max(ur),"d",cex=1.5)

dev.copy2pdf(file="potentials.pdf")

X11()


### Figure for canonical distributions

par(mfrow=c(2,2),mar=0.5*c(4, 4, 4, 2) + 0.1)

xv <- seq(-2,2,0.01)

xr <- range(xv)
ur <- c(-0.5,2)

# Constant potential
plot(xr,ur,type="n",xlab="x",ylab=expression(phi(x)),main=expression(phi[0](x)))
lines(c(-1,1),0.5*c(1,1),lwd=3)
lines(c(min(xr),-1),0*c(1,1),lwd=3)
lines(c(max(xr),1),0*c(1,1),lwd=3)
lines(c(-1,-1),ur,lty="dashed")
lines(c(1,1),ur,lty="dashed")

text(min(xr)+0.1*diff(xr),max(ur),"a",cex=1.5)


# Linear potential
D1 <- 0.5
D2 <- 0.25
phir <- c(0,1/min(D1,D2))
plot(xr,phir,type="n",xlab="x",ylab=expression(phi[1](x)),main=expression(phi[1](x)))

xp <- seq(0,max(xr),length=101)
lines(c(min(xr),0),c(0,0),lwd=3)
lines(c(0,0),phir,lty="dashed")

lines(xp,1/D1*exp(-xp/D1),lwd=3)
lines(xp,1/D2*exp(-xp/D2),lwd=1)

text(min(xr)+0.1*diff(xr),max(ur),"b",cex=1.5)

## Quadratic potential

phir <- c(0,1.5)

D2 <- 0.1
plot(xr,phir,type="n",xlab="x",ylab=expression(phi(x)),main=expression(phi[2](x)))
lines(xv,dnorm(xv,mean=0,sd=sqrt(D1)),lwd=3)
lines(xv,dnorm(xv,mean=0,sd=sqrt(D2)),lwd=1)

text(min(xr)+0.1*diff(xr),max(ur),"c",cex=1.5)


# Double-well
plot(xr,phir,type="n",xlab="x",ylab=expression(phi[4](x)),main=expression(phi[4](x)))

phi41 <- exp(-U4/D1)
phi42 <- exp(-U4/D2)

Z41 <- sum(phi41)*dx
Z42 <- sum(phi42)*dx

lines(xv, phi41/Z41,lwd=3)
lines(xv, phi42/Z42,lwd=1)

text(min(xr)+0.1*diff(xr),max(ur),"d",cex=1.5)

dev.copy2pdf(file="canonical.pdf")
