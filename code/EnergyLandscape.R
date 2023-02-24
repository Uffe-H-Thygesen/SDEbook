### Figure to draw energy landscape, for the Lyapunov chapter

Xvec <- seq(-1.25,1.25,length=101)
Vvec <- seq(-2.5,2.5,length=51)

E <- outer(Xvec,Vvec,function(x,v) 1-cos(2*pi*x)+0.5*v^2)

par(cex=1,mar=c(0,1,0,0))
pdf(width=5,height=5,file="Energy-contour.pdf")

contour(x=Xvec,y=Vvec,z=E,xlab="x",ylab="v",nlevels=5)
#        main="Contour plot of total energy E(x,v)")

points(pi*(0.5+c(0,2)),numeric(2),pch=1)
points(pi*(1.5+c(0,2)),numeric(2),pch=16)

dev.off()

pdf(width=5,height=5,file="Energy-surf.pdf")

persp(Xvec,Vvec,E,xlab="x",ylab="v",zlab="E",
      d=1e6,theta=30, phi=30, expand=0.6,
      col='lightblue', shade=0.75, ltheta=120,
      ticktype='detailed')


#        main="Surf plot of total energy E(x,v)
dev.off()

