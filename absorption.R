## Verification of the formulas for absorped Brownian motion on [0,1]

require(SDEtools)

u <- function(x) 0
D <- function(x) 0.5

dx <- 0.01
xv <- seq(0,1,dx)
xc <- c(0,0.5*(head(xv,-1)+tail(xv,-1)),1)

G <- fvade(u,D,xv,'e')

## Plot probability of having been absorbed at x=0, as function of t
i0 <- round(length(xv) / 4)
x0 <- xc[i0]

dt <- 0.25
tvec <- seq(0,3,dt)^2

P <- sapply(tvec,function(t) expm(G*t)[i0,])
PY0 <- P[1,]

## Fouerier methods
N <- 16
nv <- 1:N
cn <- 2*sin(nv*pi*x0)
lambda <- -nv^2*pi^2/2


F <- outer(nv,xc,function(n,x) sin (n*pi*x))

PY0f <- sapply(tvec,function(t)sum(-0.5*cn*(1-exp(lambda*t))*nv*pi/lambda))

plot(tvec,PY0,ylim=c(0,1))
abline(h=1-x0)
lines(tvec,PY0f)

