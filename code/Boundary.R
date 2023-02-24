## Plot of the probability of exit to the right

PlotPXtaub <- function(mu,sigma2,a,b,...)
    {
        eta <- 1- 2*mu/sigma2
        h <- function(x) (x^eta - a^eta)/( b^eta - a^eta)

        plot(h,from=a,to=b,...)
    }

sigma2 <- 1

pdf(file="Boundary.pdf",width=8,height=4)

par(mfrow=c(1,2))
    
mu <- 1

PlotPXtaub(mu,sigma2,1e-8,1,xlim=c(0,1),ylim=c(0,1),lty=1,lwd=2,main=expression(mu==1*", "*sigma^2*"="*1))
PlotPXtaub(mu,sigma2,0.001,1,add=TRUE,lty=2)
PlotPXtaub(mu,sigma2,0.01,1,add=TRUE,lty=3)
PlotPXtaub(mu,sigma2,0.1,1,add=TRUE,lty=4)

legend("bottomright",lty=c(1,2,3,4),legend=c("a=1e-8","a=0.001","a=0.01","a=0.1"),lwd=c(2,1,1,1))

mu <- 0.25

PlotPXtaub(mu,sigma2,1e-8,1,xlim=c(0,1),ylim=c(0,1),lty=1,lwd=2,main=expression(mu==0.25*", "*sigma^2*"="*1))
PlotPXtaub(mu,sigma2,0.001,1,add=TRUE,lty=2)
PlotPXtaub(mu,sigma2,0.01,1,add=TRUE,lty=3)
PlotPXtaub(mu,sigma2,0.1,1,add=TRUE,lty=4)

legend("bottomright",lty=c(1,2,3,4),legend=c("a=1e-8","a=0.001","a=0.01","a=0.1"),lwd=c(2,1,1,1))


dev.off()
