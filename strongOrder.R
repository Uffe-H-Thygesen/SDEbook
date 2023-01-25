## Compute the strong order of a scheme by comparing to the finest discretization

require(SDEtools)

strongOrder <- function(scheme,f,g,x0,T,dt,ndouble,nrepl)
{
    XT <- array(NA,c(length(x0),nrepl))
    XdT <- array(NA,c(length(x0),ndouble,nrepl))
    XeT <- array(NA,c(length(x0),ndouble,nrepl))
    XneT <- array(NA,c(ndouble,nrepl))
                  
    tsim <- seq(0,T,dt)

    for(i in 1:nrepl)
    {
        B <- rvBM(tsim,length(x0))

        sim <- scheme(f,g,tsim,x0,B)

        XT[,i] <- sim$X[length(tsim),]

        tv <- tsim
        
        for(j in 1:ndouble)
        {
            I <- seq(1,length(tv),2)
            tv <- tv[I]
            B <- matrix(B[I,],nrow=length(tv))

            simr <- scheme(f,g,tv,x0,B)
            XdT[,j,i] <- simr$X[length(tv),]
            XeT[,j,i] <- simr$X[length(tv),] - XT[,i]
            XneT[j,i] <- sum(XeT[,j,i]^2)
        }
    }

    h <- dt* 2^seq(1,ndouble)
    e <- sqrt(apply(XneT,1,mean))
    fit <- lm(log(e) ~ log(h))
    slope <- fit$coef[2]
    
    plot(h,e,log="xy",xlab="Time step",ylab="R.m.s. error")

    lines(h,exp(predict(fit)))

    legend("bottomright",legend=sprintf("%.1g power law",slope))

    return(list(XT=XT,XdT=XdT,XeT=XeT,XneT=XneT))
}

## strongOrder(euler,function(x) x,function(x) x,1,1,2^-10,4,1000)
## dev.new()
## strongOrder(euler,function(x) x,function(x) x,1,1,2^-10,4,1000)
## dev.new()

heunP <- function(...) heun(...,p=abs)

r <- 0.1
lambda <- 0.2
xi <- 1
s <- 1

f <- function(xs) c(r*xs[1],lambda*(xi-xs[2]))
g <- function(xs) diag(c(xs[1]*sqrt(xs[2]),s*xs[2]))

system.time(res <- strongOrder(heun,f,g,c(1,1),1,2^-10,5,10000))
