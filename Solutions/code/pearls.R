## Stochastic simulation of n pearls on a string

require(Matrix)
require(SDEtools)
require(MASS)

## Parameters
n <- 6
c <- 0.1
kappa <- 2
sigma <- 1

## System
K <- bandSparse(n,n,
                c(-1,0,1),
                list(rep(-kappa,n-1),rep(2*kappa,n),rep(-kappa,n-1)))

O <- sparseMatrix(i=numeric(0),j=numeric(0),dims=c(n,n))
I <- Diagonal(n,x=1)

A <- rbind(cbind(O,I),cbind(-K,-c*I))
G <- rbind(O,I)

## Initial condition for the simulation
x0 <- numeric(n)
v0 <- numeric(n)

## alternatively, sample from the stationary distribution
## v0 <- rnorm(n,sd=sigma^2/2/c)
## x0 <- mvrnorm(mu=numeric(n),Sigma=solve(K)*sigma^2/2/c)

xv0 <- c(x0,v0)

## Simulation using the transition probabilities
dt <- 1
times <- seq(0,900,dt)
tp <- dLinSDE(A,G,dt)

XV <- array(NA,c(2*n,length(times)))
XV[,1] <- xv0
for(i in 2:length(times))
    XV[,i] <- tp$eAt %*% as.numeric(XV[,i-1]) + mvrnorm(mu=numeric(2*n),Sigma=tp$St)

## Compute the empirical covariance, removing burn-in
burnin <- 20
CXV.emp <- cov(t(XV[,times>burnin]))

## Compare with the analytical
CXV.ana <- lyap(A,G %*% t(G))

## To compare the two, use for example the cosine to the angle between the two
tr <- function(A) sum(diag(A))
print(tr(CXV.emp %*% CXV.ana) / norm(CXV.ana,"F") /norm(CXV.emp,"F"))
## (1 would indicate identical matrices up to scale, etc)
print(tr(CXV.emp)/tr(CXV.ana))
## (1 would indicate identical scale)

## Plot at times
plot.times <- seq(0,900,100)
plot.ind <- times %in% plot.times

matplot(plot.times,t(XV[1:n,plot.ind]),type="l")

## Note: The simulation can also be done with the Euler-Maruyama method, or the Heun method, but this would require excessively short time steps. For example, the Euler-Maruyama method with a time step of 0.01:
sim.times <- seq(0,9,0.01)
f <- function(x) A %*% x
g <- function(x) G
print(system.time(sim.euler <- euler(f,g,sim.times,xv0)))
CXV.euler <- cov(sim.euler$X[sim.times>burnin,])
print(tr(CXV.euler)/tr(CXV.ana))
## Note that the sum-of-variances is much bigger for the simulation than
## for the analytical result. This is due to the time discretization.
## It can be improved somewhat with the Heun method (note that the noise is
## additive, so the Heun method is consistent)

sim.heun <- heun(f,g,sim.times,xv0)
CXV.heun <- cov(sim.heun$X[sim.times>burnin,])
print(tr(CXV.heun)/tr(CXV.ana))


