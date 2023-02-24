## Find the distribution of the maximum of a Brownian bridge

T <- 2
b <- 1.5

# A) Analytical

P.S.ge.x <- function(x)dnorm(2*x-b,mean=0,sd=sqrt(T))/dnorm(b,mean=0,sd=sqrt(T))

xvec <- seq(b,b+3*sqrt(T),length=101)

plot(xvec,P.S.ge.x(xvec),type="l")

# Simulation

Nsample <- 1000
dt <- 0.0001

tvec <- seq(0,T,dt)

dB <- array(rnorm(numeric(Nsample*(length(tvec)-1)),mean=0,sd=sqrt(dt)),
            c(length(tvec)-1,Nsample))

B <- apply(dB,2,function(db)c(0,cumsum(db)))

BB <- apply(B,2,function(bsim)bsim+tvec/T*(b-bsim[length(bsim)]))


# plot(c(0,T),c(0-sqrt(T),b+sqrt(T)),type="n")
# apply(BB,2,function(bb)lines(tvec,bb))

S <- apply(BB,2,max)

lines(sort(S),seq(1,0,length=Nsample),type="S")
