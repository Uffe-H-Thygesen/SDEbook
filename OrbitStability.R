### Script to illustrate orbit stability

par(pty="s")

T <- 5.5
tvec <- seq(0,T,length=101)

plot.path <- function(r,...)
  {
    omega <- 1/r^1
    x <- r*cos(omega*tvec)
    y <- r*sin(omega*tvec)
    points(x[1],y[1],pch=16,col="green")
    points(x[length(x)],y[length(y)],pch=16,col="red")
    lines(x,y,...)
    return(list(x=x,y=y))
  }

plot(c(-1,1),c(-1,1),type="n",
     main="Circular orbits",
     xlab="x",ylab="y")

plot.path(1)
plot.path(0.95)
