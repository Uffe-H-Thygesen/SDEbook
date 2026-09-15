### Illustration of random walk and its Brownian limit

set.seed(123456)

### Simulation of random walk 

## Rnadom Walk Parameters
k <- 1     ## Steps in space
h <- 1     ## Time steps
p <- 0.25  ## Probability of stepping each side

## Number of steps
N <- 10000

## Random variables for deciding if and where to step
U <- runif(N)

## Generate increments in the random walk
dX <-   (U>(1-p)) - (U<p) 

## Generate random walk
X <- c(0,cumsum(dX))

T <- (0:N)*h

Nplot1 <- 100

pdf("figRW.pdf",height=4)
par(mfrow=c(1,2))

## Plot the first part of the curve
plot(T[1:Nplot1],X[1:Nplot1],type="s",xlab="Time [ns]",ylab=expression(paste("Position [",mu,"m]")))

## plot the entire curve
plot(T,X,type="s",xlab="Time [ns]",ylab=expression(paste("Position [",mu,"m]")))
dev.off()

