## Solution of the problem of optimal control for a predator-prey system

require(SDEtools)
require(fields)
require(pracma)

graphics.off()

seed <- 123456

## Generic function for solving a continuous-time Markov chain
## Can run for a specific time, or a specific number of jumps
simulate <- function(G,i0,T,Nmax=Inf)
{
    Gt <- t(G)
    
    i <- numeric(min(Nmax,1e4))
    t <- i
    j <- 1

    i[1] <- i0
    
    while(TRUE)
    {
        lambda <- -G[i[j],i[j]]

        if(lambda==0) break
        
        t[j+1] <- t[j] + rexp(n=1,rate=lambda)

        ## Extract vector of next states as a sparse vector
        v <- Gt[,i[j],drop=FALSE]

        ##
        p <- pmax(0,v@x)
        p <- p/sum(p)

        i[j+1] <- v@i[sample(length(p),size=1,replace=TRUE,prob=p)] + 1

        if(t[j+1] > T) break
        if(j>=Nmax) break

        if( (j+1) == length(i) ) {
            i <- c(i,numeric(1e4))
            t <- c(t,numeric(1e4))
        }

        j <- j+1

    }

    i <- i[1:j]
    t <- t[1:j]
    
    return(list(t=t,i=i))
}

## Computational grid 
xgrid <- seq(-4,1,length=151)
ygrid <- seq(-4,1,length=152)

nx <- length(xgrid)-1
ny <- length(ygrid)-1

xc <- 0.5*(head(xgrid,-1) + tail(xgrid,-1))
yc <- 0.5*(head(ygrid,-1) + tail(ygrid,-1))

dx <- mean(diff(xgrid))
dy <- mean(diff(ygrid))

centers <- cell.centers(xgrid,ygrid)
xx <- pack.field(centers$x)
yy <- pack.field(centers$y)

## Stochastic Rosenzweig-MacArthur model: x=log(N), y=log(P)
## Advective flux field. Note the Ito correction for diffusivity 
ux <- function(x,y) r*(1-exp(x)/K) - c*exp(y)/(exp(x)+Nbar) - Dx
uy <- function(x,y) e*c*exp(x)/(exp(x)+Nbar) - mu - Dy

## Parameters
r <- 1
K <- 1
c <- 1 # 0.25
mu <- 0.15
Nbar <- 1.5
e <- 0.5 # 0.22
    
## Diffusivities
Dx <- 0.01
Dy <- 0.01

solver <- function(do.simulate=TRUE)
{
    set.seed(seed)


    G0 <- fvade2d(ux,uy,function(x,y)Dx,function(x,y)Dy,xgrid,ygrid)

    ## Determine stationary density (w.r.t. dn*dp) without fishing 
    pi0 <- StationaryDistribution(G0)
    pi0.field <- unpack.field(pi,nx,ny)/outer(diff(exp(ygrid)),diff(exp(xgrid)))

    ## Controls: Fisheries 1 on prey, 2 on predators. The control is the mortality
    ux1 <- function(x,y) -1
    uy1 <- function(x,y) 0

    ux2 <- function(x,y) 0
    uy2 <- function(x,y) -1

    Dnull <- function(x,y) 0

    G1 <- fvade2d(ux1,uy1,Dnull,Dnull,xgrid,ygrid)
    G2 <- fvade2d(ux2,uy2,Dnull,Dnull,xgrid,ygrid)

    ## Optimal harvest. Maximize sqrt(C1)+sqrt(C2)
    p2 <- 1.00
    Fmax <- 3

    ## Payoff
    k <- function(u) p1*sqrt(u[,1]*exp(xx))+p2*sqrt(u[,2]*exp(yy))

    ## Optimal policies. Enforce 0 harvest at the lower and left boundaries.
    ## Enforce F >= Fmax everywhere. Avoid dV/dx <= 0 or dV/dy <= 0.
    u1 <- function(G1V) pmin(Fmax*(xx>min(xx)),p1^2*exp(xx)/4/pmin(G1V,-1e-10)^2)
    u2 <- function(G2V) pmin(Fmax*(yy>min(yy)),p2^2*exp(yy)/4/pmin(G2V,-1e-10)^2)

    uopt <- list(u1,u2)

    ## Identify the optimal strategy
    print("Policy iteration ...")
    sol <- PolicyIterationSingular(G0,list(G1,G2),k,uopt,do.minimize=FALSE)

    U1 <- unpack.field(sol$u[,1],nx,ny)
    U2 <- unpack.field(sol$u[,2],nx,ny)

    Gc <- G0 + Diagonal(x=sol$u[,1]) %*% G1 + Diagonal(x=sol$u[,2]) %*% G2

    Generator2drift <- function(G,f)
    {
        df <- G %*% f
        return(unpack.field(df,nx,ny))
    }

    ## dx <- Generator2drift(Gc,xx)
    ## dy <- Generator2drift(Gc,yy)

    pi.field <- unpack.field(sol$pi,nx,ny)/outer(diff(exp(ygrid)),diff(exp(xgrid)))

    res <- list(pi0=pi0.field,pi=pi.field,U1=U1,U2=U2)
    
    if(do.simulate)
    {
        print("Simulating ...")
        ## Simulate SDE using Euler and discretized MC
        ## Use that the diffusivity is constant 
        sx <- sqrt(2*Dx)
        sy <- sqrt(2*Dy)

        f <- function(xy)
        {
            x <- min(max(xy[1],xc[1]),xc[length(xc)])
            y <- min(max(xy[2],yc[1]),yc[length(yc)])
            
            u1 <- interp2(xc,yc,U1,xp=x,yp=y)
            u2 <- interp2(xc,yc,U2,xp=x,yp=y)
            
            return(c(ux(x,y) - u1,uy(x,y) - u2))
        }

        g <- function(xy) return(diag(c(sx,sy)))

        Tsim <- 500
        dt <- 0.01 
        sim <- SDEtools::euler(f,g,seq(0,Tsim,dt),x0 =c(-3,-3))

        p <- function(x,i) pmin(pmax(x,i[1]),i[2])
        sim$u1 <-  interp2(xc,yc,U1,
                              xp=p(sim$X[,1],range(xc)),
                              yp=p(sim$X[,2],range(yc)))
        sim$u2 <-  interp2(xc,yc,U2,
                              xp=p(sim$X[,1],range(xc)),
                              yp=p(sim$X[,2],range(yc)))


        ## Discretized closed-loop generator
        ## sim <- simulate(Gc,i0=1,T=100)

        sim$Cx <- exp(sim$X[,1])*sim$u1
        sim$Cy <- exp(sim$X[,2])*sim$u2

        res <- c(res,sim=sim)
    }

    return(res)
}

my.image <- function(x,y,z,col=col,xlab,ylab,main="",xl=NULL,yl=NULL)
{
    layout(matrix(1:2,nrow=1),width=c(0.75,0.25))
    require(fields)

    xa <- seq(min(x),max(x),length=101)
    ya <- seq(min(y),max(y),length=99)
    m <- meshgrid(xa,ya)

    za <- interp.surface(list(x=x,y=y,z=z),cbind(as.numeric(m$X),as.numeric(m$Y)))
    za <- array(za,c(length(ya),length(xa)))
    
    par(mar=c(5,4,4,0)+0.1)
    plot(range(x),range(y),type="n",xlab=xlab,ylab=ylab,main=main,xaxs="i",yaxs="i")
    rangex <- range(x)
    rangey <- range(y)
    rangez <- range(z)
    scale <- function(z) (z-rangez[1])/diff(rangez)
    rasterImage((scale(za[nrow(za):1,])),rangex[1],rangey[1],rangex[2],rangey[2])

    if(!is.null(xl)) lines(xl,yl)

    par(mar=c(5,0,4,4)+0.1)
    plot(c(0,1),rangez,type="n",xlab="",xaxt="n",axes=FALSE,yaxs="i")
    axis(4)
    rasterImage(seq(1,0,length=101),0,rangez[1],1,rangez[2])
}

## Plotter
plot.res <- function(res,filename)
{
    ## Plot stationary p.d.f. without fishing
    ## Note, pdf w.r.t. the natural co-ordinates, not log-transformed

    ## Plotting region
    xmax <- 1.2
    ymax <- 0.8
    
    ## Cells to be plotted
    Ix <- xc <= log(xmax)
    Iy <- yc <= log(ymax)

    ## Interfaces to be 
    pdf(file=paste(filename,"-a.pdf",sep=""),width=4,height=4)
    with(res,
         my.image(exp(xc[Ix]),
                  exp(yc[Iy]),
                  t(pi[Iy,Ix]),
                  xlab="Prey N",ylab="Predators P"))
    dev.off()
    
    It <- res$sim.times >= (tail(res$sim.times,1)-200 ) 
    It <- It & ((1:length(It)) %% 50 ==0)
    
    pdf(file=paste(filename,"-b.pdf",sep=""),width=4,height=4)
    with(res,
         plot(sim.times[It],sqrt(sim.Cx[It]),
              type="l",ylim=c(0,sqrt(max(c(sim.Cx[It],sim.Cy[It])))),
              xlab="Time t",
              ylab=expression("Rewards "*sqrt(C^N)*","*sqrt(C^P)))
         )
    legend("topright",lty=1:2,legend=c("Prey","Predator"))
    
    with(res, lines(sim.times[It],sqrt(sim.Cy[It]),lty=2) )

    dev.off()

    pdf(file=paste(filename,"-c.pdf",sep=""),width=4,height=4)
    with(res,
         contour(exp(xc[Ix]),exp(yc[Iy]),t(U1[Iy,Ix]),
                 levels=seq(0.001,0.01,0.001),
                 xlab="Prey N",ylab="Predators P",
                 main="Prey harvest")
         )
    dev.off()

    pdf(file=paste(filename,"-d.pdf",sep=""),width=4,height=4)
    with(res,
         contour(exp(xc[Ix]),exp(yc[Iy]),t(U2[Iy,Ix]),
                 levels = seq(0.02,0.06,0.01),
                 xlab="Prey N",ylab="Predators P",
                 main="Predator harvest")

         )
    dev.off()
}

## p1s <- c(0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.5,0.3,0.4,0.5)

p1s <- 0.05

res <- list()

for(i in 1:length(p1s))
{
    p1 <- p1s[i]
    resi <- solver()
    resi <- c(resi,p1=p1)
    res[[i]] <- resi
}

save.image(file="control2D.RData")

for(i in 1:length(p1s)) #  1:length(p1s))
{
    print(filename <- paste("control2D_",format(res[[i]]$p1,dec="-"),sep=""))
    plot.res(res[[i]],filename)
}

