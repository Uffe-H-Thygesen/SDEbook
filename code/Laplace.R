## Figure to illustrate the Laplace approximation

pdf(file="Laplace.pdf",width=5,height=4)
xmin <- -2
xmax <- 4

logf <- function(x) -0.5*(x-1)^2 + 1/6*(x-1)^3 - 1/16*(x-1)^4

f <- function(x) exp(logf(x))

plot(f,from=xmin,to=xmax,lwd=2)

hatx <- optimize(logf,c(xmin,xmax),maximum=TRUE)$maximum

tiny <- 1e-3

H <- - (logf(hatx+tiny) + logf(hatx-tiny) - 2*logf(hatx))/tiny^2

plot(function(x) f(hatx) * exp(-0.5*H*(x-hatx)^2),from=xmin,to=xmax,add=TRUE,lty="dashed")

(I <- integrate(f,lower = -Inf, upper= Inf))

(Ihat <- f(hatx) * sqrt(2*pi) / sqrt(H))

legend("topright",legend=c(expression(f(x)),expression(hat(f)(x))),lty=c("solid","dashed"),lwd=c(2,1))

dev.off()
