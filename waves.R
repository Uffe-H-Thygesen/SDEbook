# Demonstration of attenuation of waves by diffusion

pdf(file="waves.pdf",width=8,height=4)
par(mfrow=c(1,2))
xx <- seq(0,1,length=1001)


D <- 1
T <- 0.01

k <- 2*pi

plot(xx,sin(k*xx),type="l",xlab="x",ylab="C")

lines(xx,exp(-D*k^2*T)*sin(k*xx),type="l",lty="dashed")

# text(0.4,0.9,"C(x,0)")
# text(0.24,0.5,"C(x,T)")

k <- 4*pi

plot(xx,sin(k*xx),type="l",xlab="x",ylab="C")
lines(xx,exp(-D*k^2*T)*sin(k*xx),type="l",lty="dashed")

# text(0.23,0.9,"C(x,0)")
# text(0.12,0.1,"C(x,T)")


dev.off()
