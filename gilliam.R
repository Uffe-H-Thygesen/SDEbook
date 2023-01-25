### Script to illustrate dynamic programming: Gilliam's rule

### State space: {1,2}, 1 being dead, 2 being alive

### Control: Foraging effort Choose between a safe investment and a risky one.

### Objective: Maximize expected lifetime foraging


### Number of time steps (years)
T <- 20

### Mortality
mort <- function(u) 1-p*(1-u)
harv <- function(u) u/(u0+u)

### Optimal effort for given fitness to next
###
### Note: This function is coded by hand to match the mort() and harv() above 
ustar <- function(v) max(0,min(1,sqrt(u0/(p*v))-u0))

### Fitness for a given effort and future fitness
Vstar <- function(u,v) harv(u) + (1-mort(u))*v

p <- 0.9
u0 <- 0.6

### Set up the DP solution
U <- V <- rep(NA,T)

U[T] <- 1
V[T] <- Vstar(1,0)

### Backward iteration in the DP equation
for(t in (T-1):1)
    {
        U[t] <- ustar(V[t+1])
        V[t] <- Vstar(U[t],V[t+1])
    }


### Plot DP solution
graphics.off()

pdf(file="DP-gilliam.pdf",width=6)
pch = 16

par(mfrow=c(2,1))
plot(1:T,U,xlab="Time",ylab="Optimal effort u",pch=pch,ylim=c(0,1))
plot(1:T,V,xlab="Time",ylab="Fitness V",pch=pch)
dev.off()

### Plot Gilliam's rule
pdf(file="Tradeoff-gilliam.pdf",width=6)
UU <- seq(0,1,0.01)
plot(mort(UU),harv(UU),type="l",lwd=3,xlim=c(0,1),
     xlab=expression("Mortality "*mu(u)),ylab=expression("Harvest "*rho(u)))
points(mort(U),harv(U))
points(c(mort(U[1])),c(harv(U[1])),pch="+",cex=3)
abline(a=0,b=harv(U[1])/mort(U[1]),lty="dashed")
text(mort(UU[1]),harv(UU[1]),"u=0",adj=-0.5)
text(mort(tail(UU,1)),harv(tail(UU,1)),"u=1",adj=c(0.5,-0.5))
dev.off()
