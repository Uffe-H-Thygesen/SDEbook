## -----------------------------------------------------------------------------
  rBM <- function(t) cumsum(rnorm(length(t),mean=0,sd=sqrt(diff(c(0,t)))))

  t <- c(0,0.5,1.5,2)
  N <- 1000
  
  B <- sapply(1:N,function(i)rBM(t))
  print(apply(B,1,mean))
  print(var(t(B)))


## -----------------------------------------------------------------------------
N <- 1000
t <- seq(0,1,0.001)
B <- sapply(1:N,function(i)rBM(t))
S <- apply(B,2,cummax)
S1 <- apply(B,2,max)
matplot(t,B[,1:3],type="l")
matplot(t,S[,1:3],type="l",add=TRUE)


## -----------------------------------------------------------------------------
hist(S1,freq=FALSE)
Spdf <- function(x) 2*dnorm(x)
plot(Spdf,add=TRUE,from=0,to=max(S1))


## -----------------------------------------------------------------------------
plot(ecdf(S1))
Scdf <- function(x) 2*pnorm(x)-1
plot(Scdf,add=TRUE,from=0,to=max(S1),col="red")


## -----------------------------------------------------------------------------
  s <- 0.5
  tau <- apply(S,2,function(x)t[sum(x<s)])
  hist(tau,freq=FALSE)
  taucdf <- function(t)2-2*pnorm(s/sqrt(t))
  taupdf <- function(t)dnorm(s/sqrt(t))*s*t^(-3/2)
  curve(taupdf,add=TRUE,from=0,to=max(tau))


## -----------------------------------------------------------------------------
plot(ecdf(tau))
curve(taucdf,from=0,to=max(tau),add=TRUE,col="red")


## -----------------------------------------------------------------------------
T <- 1

N <- 2^20
h <- T/N

Ndouble <- 8

B <- rBM(seq(0,T,h))
dB <- diff(B)

var <- array(0,c(Ndouble,2))

for(i in 1:Ndouble)
    {
        var[i,1] <- sum(abs(dB))
        var[i,2] <- sum(dB^2)
        dB <- apply(array(dB,c(2,length(dB)/2)),2,sum)
    }

hs <- h*2^(1:Ndouble)/2

par(mfrow=c(1,2))
plot(log10(hs),var[,1],xlab='log10(h)',ylab="Discretized total variation V(B)",pch=16,log="y")
lines(log10(hs),sqrt(2/pi/hs))
plot(log10(hs),var[,2],xlab='log10(h)',ylab=expression("Discretized quadratic variation [B]"[1]),ylim=c(0.9,1.1),pch=16,log="y")
abline(h=1)

