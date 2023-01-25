# Code to simulate superdiffusion, for "perspectives" chapter

# Survival function
lambda <- 1
k <- 2.5
G <- function(t) 1/(1+lambda*t)^k

rsojourn <- function(n,lambda=1,k=1.5)
  return(((1/runif(n))^(1/k) -1 )/lambda)

Nrep <- 1000
N <- 1000

U <- rsojourn(Nrep*(N-1),lambda=lambda,k=k)
U1 <- rsojourn(N,lambda=lambda,k=k-1)

plot(sort(U),(length(U):1)/length(U),type="l",log="xy",xlim=c(1e-3,1e2))
grid()

D <- rbinom(Nrep*N,size=1,prob=0.5)*2-1
X <- apply(array(D*U,c(Nrep,N)),1,sum)

