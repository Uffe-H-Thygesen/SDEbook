require(SDEtools)

tv <- seq(0,10,1e-3)
B <- rBM(tv)
G <- cos(B) + sin(tv)
M <- itointegral(G,B)
X <- M - 0.5*QuadraticVariation(M)
Y <- exp(X)

plot(tv,Y,type="l")
lines(tv,1+itointegral(Y*G,B),col=2)

