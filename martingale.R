s <- 1
t <- 3

h <- 1e-4

svec <- seq(0,s,h)
dB <- rnorm(length(svec)-1,sd=sqrt(h))
B <- c(0,cumsum(dB))

X <- B^2-svec

tvec <- seq(s,t,h)

npath <- 1000
dB2 <- array(rnorm(npath*(length(tvec)-1),sd=sqrt(h)),c(length(tvec)-1,npath))

B2 <- apply(dB2,2,function(dB)c(0,cumsum(dB)) + B[length(B)])


X2 <- B2^2-tvec

ylim <- range(c(X,as.numeric(X2)))

plot(svec,X,type="l",lwd=2,xlim=c(0,t),ylim=ylim)

apply(X2,2,function(X)lines(tvec,X,col="grey"))

xx <- seq(ylim[1],ylim[2],length=100)
phi <- dchisq(xx+t,df=1)

lines(t-phi,xx)
     


X11()
hist(X2[length(tvec),],freq=FALSE)
