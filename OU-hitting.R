## Hitting time for the OU process

## Find the hitting time numerical integration

require(deSolve)
f <- function(t,x,p) list(sqrt(2*pi)*exp(t^2/2)*(pnorm(t)-0.5))

tv <- seq(0,10,0.01)

res <- ode(0,tv,f)

Etau <- res[,2]

pdf("OU-hitting.pdf",width=7,height=3.5)
par(mfrow=c(1,2),mar=c(4,5,0.2,0.2))
plot(tv[tv<=2],Etau[tv<=2],type="l",xlab="Distance l",ylab=expression(E*tau[l]))

# lines(tv,sqrt(pi/2)/tv*exp(tv^2/2),lty="dotted")
## plot(tv,Etau/(exp(tv^2/2)/tv*sqrt(pi/2)),type="l",xlab="Distance l",
##      ylab=expression(E*tau[l] / (sqrt(pi/2)*exp(l^2/2)/l) ))
## abline(h=1,lty="dashed")

plot(tv[tv>=2],Etau[tv>=2],type="l",log="y",xlab="Distance l",
     ylab=expression(E*tau[l]))
# lines(tv,sqrt(pi/2)/tv*exp(tv^2/2),lty="dotted")
dev.off()
