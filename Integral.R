### Sketch of integral

par(cex=1.25,lab=c(3,3,0))

X <- function(w) 1.5+ 3*dnorm(10*w-5) + tanh(20*w-10)

wgrid <- c(0,0.38,0.55,0.62,1)
xgrid <- numeric(length(wgrid)-1)

for(i in 2:length(wgrid))
  xgrid[i-1] <- min(X(seq(wgrid[i-1],wgrid[i],length=100))) - 0.1

maxX <- max(X(seq(0,1,length=100)))

plot(X,from=0,to=1,lwd=3,
     ylim =c(0,maxX),
     xlab=expression(omega),ylab="x")

px <- rep(wgrid,rep(2,length(wgrid)))
py <- c(0,rep(xgrid,rep(2,length(xgrid))),0)

polygon(px,py,col="grey")


text(0.9,X(0.9),expression(X(omega)),pos=3)
text(0.47,approx(wgrid,c(xgrid,0),0.5,method="constant")$y,
     expression(X[s](omega)),pos=3)

text(0.8,0.5*X(0.8),expression(E(X[s])))

dev.copy2pdf(file="Integral.pdf")
