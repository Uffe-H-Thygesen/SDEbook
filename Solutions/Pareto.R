r <- 0.25
sigma <- 1

if(r - sigma^2 /2 > 0) stop("Wrong parameters!")

N <- 1000
T <- 10/(sigma^2/2 - r)

tv <- seq(0,T,length=1001)

require(SDEtools)

B <- rvBM(tv,N)

X <- exp((r-sigma^2/2)*tv+sigma*B)

S <- apply(X,2,max)

plot(sort(S),(N:1)/N,log="xy")


lambda <- 1-2*r/sigma^2
Ss <- sort(S)
Gs <- Ss^(-lambda)

lines(Ss,Gs)
