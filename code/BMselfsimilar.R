# Figure to illustrate self-similarity of Brownian motion

require(SDEtools)

set.seed(123)

# No. sample paths
N <- 3

# Time
TimeScales <- c(0,0.1,10,1000)
Resolution <- TimeScales/1000

tvec <- 0

## Generate a vector of time points which is denser in the beginning
## so that for each TimeScale, there are 1000 points
for(i in 2:length(TimeScales))
  {
    tvec <- c(tvec,seq(TimeScales[i-1],TimeScales[i],Resolution[i])[-1])
  }

## Generate all sample paths
B <- rvBM(tvec,N)

## As range on the y-axis, choose two times the standard deviation
Bscales <- 2*sqrt(TimeScales[-1])

pdf("BMselfsimilar.pdf",width=8,height=4)

## Create a row of plots
par(mfrow=c(1,length(TimeScales)-1)) 

for(k in 1:(length(TimeScales)-1))
{
    ## Create plotting window
    plot(c(0,TimeScales[k+1]),Bscales[k]*c(1,-1),
       xaxs = "i", yaxs = "i",
       type="n",xlab="Time t", ylab=expression(B[t](omega)))

    ## Add the sample paths
    for(i in 1:N) lines(tvec,B[,i],lty=i,lwd=2)

    ## Add confidence limits
    for(i in c(-1,1)) lines(tvec,i*sqrt(tvec),lty="dashed")
}  

dev.off()

