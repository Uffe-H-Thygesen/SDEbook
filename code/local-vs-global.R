require(SDEtools)

pdf(file="local-vs-global.pdf",width=5,height=5)
T <- 3
h <- 1

set.seed(123456)

dt <- 0.01
tv <- seq(0,T,dt)
B <- rBM(tv)

r <- 1
sigma <- 0.1
x <- 1

X <- x*exp( (r-0.5*sigma^2)*tv + sigma*B)

plot(tv,X,type="l",lwd=3,xlab="t",ylab=expression(X[t]))

tc <- seq(0,T,h)
Bc <- approx(tv,B,tc)$y

f <- function(x) r*x
g <- function(x) sigma*x

Xc <- euler(f,g,tc,x,Bc)$X
lines(tc,Xc,lwd=3,lty="dashed")
points(tc,Xc,pch=16)

tr <- tv
Xr <- X


for(i in 2:length(tc))
{
    Xt <- approx(tr,Xr,tc[i])$y
    XT <- approx(tr,Xr,T)$y
    
    tr <- seq(tc[i],T,dt)
    Br <- approx(tv,B,tr)$y
    Xr <- Xc[i]*exp( (r-0.5*sigma^2)*(tr-tc[i]) + sigma*(Br-Br[1]))
    lines(tr,Xr)
    arrows(T,XT,T,tail(Xr,1),length=0.1,lwd=2,col="grey")
    arrows(tc[i],Xt,tc[i],Xc[i],length=0.1)
    points(rep(T,2),c(XT,tail(Xr,1)),pch=16,cex=0.5)
}

dev.off()
