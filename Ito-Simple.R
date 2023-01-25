#  Three panels, or three curves in one panel: a) A sample path of Brownian motion. b) A sample path of an elementary function. c) The corresponding sample path of $\int_S^T f_t ~dB_t$, as function of $T$.

#  Perhaps force the realization so that e.g. B is negative at time S, f is positive, while the integral is around 0.


T <- 1
h <- 1e-3

t <- seq(0,T,h)

dB <- rnorm(length(t)-1,sd=sqrt(h))
B <- c(0,cumsum(dB))

G <- 1 + 0.25*sin(2*pi*t) + 0.25*B


ts <- t[seq(1,length(t),100)]
Bs <- B[seq(1,length(t),100)]
Gs <- F[seq(1,length(t),100)]

iGdB <- c(0,cumsum(Gs[-length(Gs)]*diff(Bs)))

par(mfrow=c(2,1))
plot(t,B,type="l")
lines(ts,iGdB)

plot(t,F,type="l")




