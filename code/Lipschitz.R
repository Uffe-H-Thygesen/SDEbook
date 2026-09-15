# Graph to illustrate Lipschitz continuity

pdf(file="Lipschitz.pdf",width=7,height=3.5)

## Lipschitz constant
K <- 4

## Range of x axis
x0 <- -5
x2 <- 5

## The Lipschitz continuousfFunction
f <- function(x) x+sin(x)+cos(2*x)

plot(f,from=x0,to=x2)

## The test point
x1 <- -1
f1 <- f(x1)

## Bound on increments in f
df0 <- (x1-x0)*K
df2 <- (x2-x1)*K

## Plot the "permitted regions"
polygon(c(x1,x2,x2),f1+c(0,1,-1)*df2,col="lightgray",border=NA)
polygon(c(x1,x0,x0),f1+c(0,1,-1)*df0,col="lightgray",border=NA)

## Plot the "forbidded zones"
polygon(c(x1,x0,x2),f1+c(0,df0,df2),col="darkgray",border=NA)
polygon(c(x1,x0,x2),f1-c(0,df0,df2),col="darkgray",border=NA)

## Add the graph again
plot(f,from=x0,to=x2,add=TRUE,lwd=4)

dev.off()
