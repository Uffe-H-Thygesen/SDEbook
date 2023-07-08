## Hitting time for the OU process
##
## dX = - X*dt  + sqrt(2)*dB

## Find the hitting time numerical integration, using an ODE solver
require(deSolve)
f <- function(t,x,p) list(sqrt(2*pi)*exp(t^2/2)*(pnorm(t)-0.5))

tv <- seq(0,10,0.01)

res <- ode(0,tv,f)

Etau <- res[,2]

## Make plot with two panels, differing only in axis
pdf("OU-hitting.pdf",width=7,height=3.5)
par(mfrow=c(1,2),mar=c(4,5,0.2,0.2))
plot(tv[tv<=2],Etau[tv<=2],type="l",xlab="Distance l",ylab=expression(E*tau[l]))

plot(tv[tv>=2],Etau[tv>=2],type="l",log="y",xlab="Distance l",
     ylab=expression(E*tau[l]))

## For comparison, we could add the approximation
## points(tv,sqrt(pi/2)/tv*exp(tv^2/2))

dev.off()
