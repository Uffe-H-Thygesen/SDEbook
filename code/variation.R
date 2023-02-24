require(latex2exp)

pdf(file="variation.pdf",width=7,height=4)

T <- 1

N <- 2^20
h <- T/N

Ndouble <- 8

dB <- rnorm(N,sd=sqrt(h))


var <- array(0,c(Ndouble,2))

for(i in 1:Ndouble)
    {
        var[i,1] <- sum(abs(dB))
        var[i,2] <- sum(dB^2)
        dB <- apply(array(dB,c(2,length(dB)/2)),2,sum)
    }

hs <- h*2^(1:Ndouble)/2

par(mfrow=c(1,2))
plot(log10(hs),var[,1],xlab=TeX("$\\log_{10}(\\Delta t)$"),ylab="Discretized total variation V(B)",pch=16,log="y")
plot(log10(hs),var[,2],xlab=expression(log[10](Delta*t)),ylab=expression("Discretized quadratic variation [B]"[1]),ylim=c(0.9,1.1),pch=16,log="y")

dev.off()
