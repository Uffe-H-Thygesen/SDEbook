require(SDEtools)

pdf(file="dw1d.pdf")

par(bg="#140E0C",fg="#FFFFFF")

f <- function(x) x-x^3
g <- function(x) sigma
D <- function(x) 0.5*sigma^2

xmax <- 1.5

sigma <- 0.1
T <- 5.9
dt <- 0.01

tv <- seq(0,T,dt)

set.seed(123456)

x0 <- 0.00

sim <- heun(f,g,tv,x0)

plot(sim$times,sim$X,type="l",ylim=c(-1,1)*xmax,lwd=1,
     xlab="",ylab="",xaxt="n",yaxt="n")


xi <- seq(-xmax,xmax,length=201)

h <- 1
tt <- seq(0.00001,T,h)

G <- fvade(f,D,xi,'a')
P <- expm(G*h)

s0 <- 0.1

phi0 <- diff(pnorm(xi,x0,s0))

PHI <- array(NA,c(length(tt),length(phi0)))
PHI[1,] <- phi0
for(i in 2:nrow(PHI)) PHI[i,] <- as.numeric(PHI[i-1,] %*% P)

for(i in 1:nrow(PHI))
{
    scale <- 1/max(PHI[i,])
    lines(tt[i]+scale*PHI[i,],xi[-1],lwd=3)
}

dev.off()
