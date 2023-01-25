pdf(file="law-iter-log.pdf",width=7,height=4)

t <- c(0,10^seq(0,250,length=10000))

dt <- diff(t)

dB <- rnorm(length(dt),sd=sqrt(dt))

B <- c(0,cumsum(dB))

plot(log(t),B / 2/sqrt(t),type="l",ylim=c(-3,3),
     xlab=expression(log*" t"),ylab=expression(B[t]/sqrt(t)))

lines(t,+2*sqrt(log(log(t))),lwd=2,col="grey")
lines(t,-2*sqrt(log(log(t))),lwd=2,col="grey")

dev.off()

# dev.copy2pdf(file="law-iter-log.pdf")
