# Figure to illustrate self-similarity of Brownian motion

# require(grid)
# require(lattice)

pdf("BMselfsimilar.pdf",width=8,height=4)

set.seed(123)

# No. sample paths
N <- 3

# Time
TimeScales <- c(0,0.1,10,1000)
Resolution <- TimeScales/1000

tvec <- 0

for(i in 2:length(TimeScales))
  {
    tvec <- c(tvec,seq(TimeScales[i-1],TimeScales[i],Resolution[i])[-1])
  }

dtvec <- diff(tvec)

Ndt <- length(dtvec)

# Simulate BM
dB <- array(0,c(Ndt,N))
for(j in 1:N)
  dB[,j] <- sqrt(dtvec)*rnorm(Ndt)

B <- apply(dB,2,function(db) c(0,cumsum(db)))

rB <- range(B)

Bscales <- 2*sqrt(TimeScales[-1])

par(mfrow=c(1,length(TimeScales)-1)) # ,cex=2,lab=c(1,2,5))
for(k in 1:(length(TimeScales)-1))
{
  plot(TimeScales[k+c(0,1)],Bscales[k]*c(1,-1),
       xlim=c(0,TimeScales[k+1]),
       xaxs = "i", yaxs = "i",
       type="n",xlab="Time t", ylab=expression(B[t](omega)))

  for(i in 1:N) lines(tvec,B[,i],lty=i,lwd=2)
  for(i in c(-1,1)) lines(tvec,i*sqrt(tvec),lty="dashed")
}  

dev.off()

