require(SDEtools)


t <- c(0,10^seq(0,250,length=10000))
B <- rBM(t)

pdf(file="law-iter-log.pdf",width=7,height=4)

plot(log(t),B / 2/sqrt(t),type="l",ylim=c(-3,3),
     xlab=expression(log*" t"),ylab=expression(B[t]/sqrt(t)))

lines(t,+2*sqrt(log(log(t))),lwd=2,col="grey")
lines(t,-2*sqrt(log(log(t))),lwd=2,col="grey")

dev.off()
