% Simulation of Bernouilli process

N <- 20
p <- 0.25
X <- rbinom(N,1,p)

plot(X,pch=16,yaxp=c(0,1,1),xaxp=c(1,N,N-1),xlab="Time i",ylab="X")

dev.copy2pdf(file="simBernouilli.pdf",width=8,height=3)
