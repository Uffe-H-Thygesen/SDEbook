## Simulations

set.seed(20)

sigma <- 1/sqrt(2) # Noise in population dynamics

T <- 100
dt <- 0.01

require(SDEtools)

tvec <- seq(0,T,dt)
BM <- rBM(tvec)

p <- function(x) max(0,x)

Xol <- euler(function(x)(x*(1-x)-1/4)*(x>0),function(x) sigma*x*(x>0),tvec,1,BM,p=p)
print(mean(sqrt(1/4*(Xol$X[,1]>0))))

Xce <- euler(function(x)(x*(1-x)-1/2*x)*(x>0),function(x) sigma*x*(x>0),tvec,1,BM,p=p)
print(mean(sqrt(1/2*Xce$X)))
print(mean(Xce$X))

Xoe <- euler(function(x)(x*(1-x)-x^2)*(x>0),function(x) sigma*x*(x>0),tvec,1,BM,p=p)
print(mean(Xoe$X))

pdf(file="Fisheries-sim.pdf",width=6,height=4.5)
par(mfrow=c(2,1),mar=c(1, 5, 3, 1)+0.1)
lty <- c(2,1,1)
lwd <- c(1,1,3)

matplot(tvec,cbind(Xol$X,Xce$X,Xoe$X),type="l",
        xlim=c(0,20),lty=lty,lwd=lwd,col="black",xaxt="n",
        xlab="",ylab=expression("Biomass "*X[t]))
legend("topright",col="black",
       lty=lty,lwd=lwd,
       legend=c("Constant harvest","Constant effort","Optimal policy"))

par(mar=c(4, 5, 0, 1)+0.1)
matplot(tvec,cbind(1/2*(Xol$X>0),sqrt(1/2*Xce$X),Xoe$X),type="l",
        xlim=c(0,20),lty=lty,lwd=lwd,,col="black",
        xlab="Time t",ylab=expression("Profit rate "*sqrt(U[t])))
dev.off()

