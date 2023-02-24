### Sample path of Ito integral of an elementary function

require(SDEtools)

set.seed(1)

plot.elementary <- function(tv,gv,...)
    {
        plot(range(tv),range(gv),type="n",...)
        for(i in 1:(length(tv)-1))
            lines(tv[i:(i+1)],rep(gv[i],2))
        points(tail(tv,-1),gv)
        points(head(tv,-1),gv,,pch=16)
    }

tv <- c(0,0.5,2,3,4,5)
gv <- c(1.5,1,0,-1,1.5)

tsim <- seq(0,tail(tv,1),0.01)

B <- rBM(tsim)

## Really, fix the endpoint of the BM

BT <- -1

B <- B + (BT-tail(B,1))*tsim/tail(tsim,1)

S <- 0.75

i <- tsim>=S
tt <- tsim[i]
g <- approx(tv,c(gv,tail(gv,1)),tt,method="constant",f=0)

I <- stochint(g$y,B[i])

pdf(file="ito-integral-elementary.pdf")

plot.elementary(tv,gv,ylim=range(c(0,gv,I,B)),xlab="T",ylab="",cex=4)
lines(tsim,B,lty="dashed")
lines(tt,I$l,lwd=3)
legend("bottomleft",c(expression(G[T]),expression(B[T]),expression(I[T])),lty=c("solid","dashed","solid"),lwd=c(1,1,3))

dev.off()
