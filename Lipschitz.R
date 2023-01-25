# Graph to illustrate Lipschitz

pdf(file="Lipschitz.pdf",width=7,height=3.5)

K <- 4

x0 <- -5
x2 <- 5

f <- function(x) x+sin(x)+cos(2*x)

plot(f,from=x0,to=x2)

x1 <- -1
f1 <- f(x1)

df0 <- (x1-x0)*K
df2 <- (x2-x1)*K

polygon(c(x1,x2,x2),f1+c(0,1,-1)*df2,col="lightgray",border=NA)
polygon(c(x1,x0,x0),f1+c(0,1,-1)*df0,col="lightgray",border=NA)

polygon(c(x1,x0,x2),f1+c(0,df0,df2),col="darkgray",border=NA)
polygon(c(x1,x0,x2),f1-c(0,df0,df2),col="darkgray",border=NA)

plot(f,from=x0,to=x2,add=TRUE,lwd=4)
dev.off()
