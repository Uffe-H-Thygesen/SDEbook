## Plot of discretized total and quadratic variation of Brownian motion
require(SDEtools)
set.seed(123)

T <- 1        ## Terminal time
N <- 2^20     ## Number of time steps, finest resolution
h <- T/N      ## Finest time step

## Number of doublings of the time step
Ndouble <- 8

## Generate the sample path with fine resolution
B <- rBM(seq(0,T,h))

## Array for storing the results; total and quadratic variation
var<- array(NA,c(Ndouble,2))

for(i in 1:Ndouble)
    {
        ## Compute total and quadratic variation
        dB <- diff(B)
        var[i,1] <- sum(abs(dB))
        var[i,2] <- sum(dB^2)

        ## Subsample 
        B <- B[seq(1,length(B),2)]
    }

## Vector of time steps
hs <- h*2^(1:Ndouble)/2

require(latex2exp)

pdf(file="variation.pdf",width=7,height=4)

par(mfrow=c(1,2))
plot(log10(hs),var[,1],
     xlab=TeX("$\\log_{10}(\\Delta t)$"),
     ylab="Discretized total variation V(B)",pch=16,log="y")
plot(log10(hs),var[,2],
     xlab=expression(log[10](Delta*t)),
     ylab=expression("Discretized quadratic variation [B]"[1]),
     ylim=c(0.9,1.1),pch=16,log="y")

dev.off()
