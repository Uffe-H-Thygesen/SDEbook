## Produces figure swim.pdf for the optimal control problem "swim right or left"

H <- 2
V <- function(x) log(cosh(H)/cosh(x))
u <- tanh

pdf(file="swim.pdf",width=6,height=4)

par(mfrow=c(1,2))
plot(V,from=-H,to=H)
plot(u,from=-H,to=H)

dev.off()
