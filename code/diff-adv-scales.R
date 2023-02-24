# Figure illustrating the square root scaling

D <- 1
u <- 1
tt <- seq(0,4,length=1001)
xx <- sqrt(2*D*tt)
yy <- u * tt
par(cex=2)
plot(tt,xx,type="l",xlab="Time t",ylab=expression("Length scales"),lwd=2)
lines(tt,yy,type="l",lwd=2,lty="dashed")
text(3.5,2.2,"Diffusive")
text(1.5,0.5,"Advective")

dev.copy2pdf(file="diff-adv-scales.pdf")
