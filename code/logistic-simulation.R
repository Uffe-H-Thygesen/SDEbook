require(SDEtools)

r <- 1
K <- 1

x0 <- 0.01

sigmas <- rev(c(0,0.1))

tv <- seq(0,25,0.01)
B <- rBM(tv)

f <- function(x) r*x*(1-x/K)


sim <- lapply(sigmas,function(s) euler(f,function(x) s*x,tv,x0,B))

xmax <- max(sapply(1:length(sigmas),function(i)sim[[i]]$X))


pdf(file="logistic-simulation.pdf",width=6,height=5)
plot(range(tv),c(0,xmax),type="n",xlab="Time",ylab="X")
sapply(sim,function(s)lines(tv,s$X))
dev.off()
