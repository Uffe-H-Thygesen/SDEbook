## Plot of a sample path from the stochastic logistic growth model,
## including the sensitivity

require(SDEtools)

set.seed(123456)

## The stochastic logistic growth model
fX <- function(x) r*x*(1-x/K)
gX <- function(x) sigma*x

## Sensitivity
fS <- function(x,s) (r-2*x/K)*s
gS <- function(x,s) sigma*s

## Extended system, containing both the original system
## and the sensitivity
fE <- function(xs) c(fX(xs[1]),fS(xs[1],xs[2]))
gE <- function(xs) sigma * xs

## Model parameters
r <- 1
K <- 1
sigma <- 0.25

## Simulation control
dt <- 0.01
T <- 15

tv <- seq(0,T,dt)

## Initial condition (for the state and for the sensitivity)
x0 <- 0.01
xs0 <- c(x0,1)

## Solve the coupled system
sol <- euler(fE,gE,tv,xs0)

pdf(file="lyap-logistic.pdf", width=8,height=4)
par(mfrow=c(1,2))
plot(tv,sol$X[,1],type="l",xlab="Time",ylab=expression(X[t]),ylim=c(0,max(sol$X[,1])))
plot(tv,sol$X[,2],type="l",xlab="Time",ylab=expression(S[t]),log="y",yaxt="none")
ax <- 10^(-2:1)
axis(2,at=ax,labels=sapply(ax,format))
dev.off()
