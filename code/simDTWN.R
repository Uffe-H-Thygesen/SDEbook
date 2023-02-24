### Simulation of Discrete Time White Noise

graphics.off()
N <- 20

X <- rnorm(N)
plot(X,pch=16,cex=2,xlab="Time i",ylab=expression(X[i]),ylim=c(-1,1)*max(abs(X)))
Xe <- c(0,X,0)

lines(rep(0:(N+1),rep(3,N+2)),rbind(rep(0,N+2),Xe,rep(0,N+2)))

dev.copy2pdf(file="simDTWN.pdf",width=8,height=5)

### Generate cumulative

S <- c(0,cumsum(X))

X11()
plot(0:N,S,pch=16,cex=2,xlab="Time t",ylab=expression(bar("B")[t]),ylim=c(-1,1)*max(abs(S)))
lines(0:N,S)
dev.copy2pdf(file="simCSWN.pdf",width=8,height=5)

