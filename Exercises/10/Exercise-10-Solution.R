## Start from a clean slate
rm(list=ls())
graphics.off()

#require(MASS)
#require(Matrix)
#require(fields)
#source("fvade.R")

require(SDEtools)


lambda <- 1
xi <- 1
gamma <- 1

# Define model. Note abs to handle negative x.
f = function(x) lambda*(xi-x);
g = function(x) gamma*sqrt(x);

# Advection-diffusion form
D = function(x) 0.5*gamma^2*x;
u = function(x) f(x) - 0.5*gamma^2;

## Generator

## Define grid
xii = seq(0,5,0.1)                                # Cell Interfaces
dx <- diff(xii)
xc = 0.5*(tail(xii,-1) + head(xii,-1))   # Cell centers

G = fvade(u,D,xii,'r');

tsample <- 0.1
P <- as.matrix(expm(G*tsample))

## Convert probabilities to densities
p2d <- function(P) apply(P,1,function(p)p/dx)

contour(xc,xc,z=P)

obs <- read.table("hmm-obs.txt",header=TRUE)

plot(Y ~t, data=obs)

vsample <- 0.5
dl <- function(x,y) dpois(x=y,lambda=vsample*x)
dltab = outer(xc,0:max(obs$Y),dl);

matplot(xc,dltab,type="l")

## State likelihood
ltab <- dltab[,obs$Y+1]


## Filter
hmmfilter <- function(G)
{
    phi <- Matrix(array(0,c(length(obs$t),length(xc))))
    psi <- phi

    ## Initialize with the stationary distribution
    mu <- StationaryDistribution(G) 
    mu <- mu/sum(mu)

    ## Compute transition probabilities
    P <- expm(G*tsample)
    const <- numeric(length(obs$t))

    phi[1,] <- mu

    ## Include the first data update
    psi[1,] <- phi[1,] * ltab[,1]
    const[1] <- sum(psi[,1] )
    psi[,1] <- psi[,1] / const[1]

    ## Main time loop over the remaining time steps
    for(i in 2:length(obs$t))
    {
        phi[i,] = psi[i-1,] %*% P    # Time update
        psi[i,] = phi[i,] * ltab[,i] # Data update
        const[i] <- sum(psi[i,])     # Normalization
        psi[i,] = psi[i,] / const[i]
    }

    return(list(c=const,phi=phi,psi=psi,loglik=sum(log(const))))
}

## Run the filter with default values for parameters    
est <- hmmfilter(P)

## Show posterior probabilities, and posterior mean 
image(obs$t,xc,t(p2d(est$psi)),zlim=c(0,2))
lines(obs$t,apply(est$psi,1,function(p)sum(p*xc)))

## Add true state 
sim <- read.table("hmm-states.txt",header=TRUE)
lines(sim$t,sim$X,lwd=2)


## Likelihood estimation of xi
xis <- seq(0,1.5,0.05)   # Try these values of xi

loglik <- function(xi)    # Likelihood function of one value of xi 
{
    f <- function(x) lambda*(xi-x);
    D <- function(x) 0.5*gamma^2*x;
    u <- function(x) f(x) - 0.5*gamma^2;


    G = fvade(u,D,xii,'r');

    return(hmmfilter(G)$loglik)
}

## Tabulate likelihood function
print(system.time(ls <- sapply(xis,loglik)))

plot(xis,ls,type="b")
grid()
abline(v=xi) ## Include the true value used in the simulation 
abline(h=max(ls)-0.5*qchisq(0.95,df=1)) ## Include confidence interval

print(system.time(opt <- optim(par=1,fn=function(xi)-loglik(xi),lower=0,upper=2)))
