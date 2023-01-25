### Sketch for the situation where a diffusion is not
### a strong Markov process

### Parameterization of curve
f <- function(b) b - 20*b*dnorm(b)
g <- function(b) dnorm(b)

### Parameter region
bvec <- seq(-10,10,0.01)

### Plot curve
plot(f(bvec),g(bvec),type="l",
     xlab="x",ylab="y")

### Add initial point
points(f(0),g(0),pch=16,cex=2)
text(f(0),g(0)-0.01,expression(X[0]),pos=1)

### Find and plot stopping point
b1 <- bvec[sum(cumprod(f(bvec)<0))]
points(f(b1),g(b1),pch=16,cex=2)
text(f(b1),g(b1)-0.01,expression(X[tau]),pos=1)

### Export
dev.copy2pdf(file="StrongMarkov.pdf")

