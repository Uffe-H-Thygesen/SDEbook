require(fields)
require(SDEtools)
require(MASS)

graphics.off()
rm(list=ls())

### Optimal fisheries management

### Parameters
K <- 1 # Carrying capacity
r <- 1 # Growth rate 
sigma <- 1/sqrt(2) # Noise in population dynamics

### Upper bound on fisheries mortality
Umax <- 5

### Discretization of state space
Xmax <- 4
dx <- 0.05
xi <- seq(0,Xmax,dx)
xc <- xi[-1] - 0.5*diff(xi)


### Functions entering into the model

### Uncontrolled system:
f <- function(x) r*x*(1-x/K)

D <- function(x) 1/2*sigma^2*x^2
dD <- function(x) sigma^2*x

u <- function(x) f(x) - dD(x)
G <- fvade(u,D,xi,'r')

### Effect of the fishing: The "generator" d/dx
ddxl <- -fvade(function(x)-1,function(x)0,xi,'r')

ddxl[1,] <- ddxl[2,]

### Simulation control
dt <- 1
T <- 10

tvec <- seq(0,T,dt)

V <- array(NA,c(length(xc),length(tvec)))
U <- V

V[,length(tvec)] <- 0
U[,length(tvec)] <- 0

tiny <- 1e-3

### Dynamic programming solution 

V <- array(NA,c(length(xc),length(tvec)))
V[,length(tvec)] <- 0
U <- V

require(expm)

## For reference: Zero and identity array
O <- array(0,c(length(xc),length(xc)))
I <- diag(rep(1,length(xc)))

for(i in (length(tvec)-1):1)
    {
        dVdx <- ddxl %*% V[,i+1]

        ## Identify optimal control:
        ## max(sqrt(x*u) - dvdx*x*u) has stationary point
        ## 1/(2*sqrt(u)) = dvdx*sqrt(x)
        ## or:  
        ustar <- 1/4/pmax(dVdx^2*xc,tiny)
        ustar <- pmax(0,pmin(ustar,Umax))
        ustar[1] <- 0

        U[,i] <- ustar

        ## Closed loop generator
        Gcl <- G - as.numeric(ustar*xc) * ddxl
        
        ## Euler stepping of value function
#        LV <- Gcl %*% V[,i+1]
#        V[,i] <- V[,i+1] + LV*dt +sqrt(ustar)*dt
        
        ## Expm stepping of value function
        GG <- rbind(cbind(Gcl,diag(sqrt(ustar*xc))),cbind(O,G))
        V[,i] <- (expm(GG*dt) %*% c(V[,i+1],rep(1,length(xc))))[1:length(xc)]
    }

U0 <- U[,1]
V0 <- V[,1]
V0 <- V0 - min(V0) 

pdf(file="fisheries-management.pdf")

par(mfrow=c(2,1),mar=c(5, 5, 1, 1)+0.1)
plot(xc,V0,xlab="Biomass x",ylab=expression("Value"*V[0](x)),type="l",lwd=1) # ,log="xy")
lines(xc,0.5*(log(xc)-log(xc[1])),lwd=3)

plot(xc,U0*xc,type="l",col="black",xlab="Biomass x",lwd=1,ylab=expression("Optimal harvest rate "*mu^"*"*(x)),log="")
lines(xc,xc^2,lwd=3)

dev.off()

