## Monte Carlo to verify the hitting time of BM with drift

require(SDEtools)

N <- 1000
u <- -0.4

T <- 10/(abs(u)+0.1)
dt <- 1e-3

times <- seq(0,T,dt)

XX <- rvBM(times,n=N,u=u)

findtau <- function(y)  match(TRUE,y>1)
tau <- times[apply(XX,2,findtau)]

ftau <- function(t)t^{-3/2}*dnorm(1/sqrt(t))*exp(u-0.5*u^2*t)

PtauFinite <- integrate(ftau,lower=0,upper=Inf)$value
Etau <- integrate(function(t)t*ftau(t),lower=0,upper=Inf)$value

tau[is.na(tau)] <- T

tab <- array(c(PtauFinite,Etau,if(u>0) 1 else exp(2*u),if(u>0) 1/u else NA,mean(tau<T),mean(tau)),c(2,3))
rownames(tab) <- c("PtauFinite","Etau")
colnames(tab) <- c("Num","Theory","Sim")
print(tab)


hist(tau,freq=FALSE)
plot(ftau,from=0,to=max(tau,na.rm=TRUE),add=TRUE)
