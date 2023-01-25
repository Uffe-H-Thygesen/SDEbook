### Plot to explain that E X = int(G,dx)

pdf(file="figSurvive.pdf")
G <- function(x) exp(-x)

xmax <- 3
xvec <- seq(0,xmax,0.01)
Gvec <- G(xvec)

par(cex=1.7)
plot(xvec,Gvec,type="l",xlab=expression(x*" or "*G^-1*(omega)),
     ylab=expression(G(x)*" or "*omega),lwd=3,ylim=c(0,1))

polygon(c(xvec,xmax,0),c(Gvec,0,0),border=NA,col="grey")

dev.off()

