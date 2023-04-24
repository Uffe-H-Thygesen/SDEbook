### Figure to illustrate the hitting time for B.M.

### What level to hit?
x <- 1

## The theoretical p.d.f. of the hitting time
density <- function(t) dnorm(-x/sqrt(t))*x*t^(-3/2)


pdf(file="BM-HittingTime.pdf",width=16,height=10)
par(mfrow=c(1,2),cex=2)

## plot the initial part of the curve
plot(density,from=1e-4,to=2,lwd=3,xlab="Time t",ylab=expression("Pdf "*f[tau](t)))
grid()

## Plot the tail of the curve
plot(density,from=1,to=1000,
     log="xy",lwd=3,xlab="Time t",ylab=expression("Pdf "*f[tau](t)))
tsup <- c(10,1000)
psup <- tsup^-1.5
lines(tsup,psup,lwd=8)
text(50,50^-1.5,"Slope = -3/2",pos=4)
grid()

dev.off()

