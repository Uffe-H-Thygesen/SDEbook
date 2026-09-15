## -----------------------------------------------------------------------------
require(SDEtools)

lambda <- 1
xi <- 1
gamma <- 1

# Define model. 
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
nc <- length(xc)

G = fvade(u,D,xii,'r');

tsample <- 0.1
P <- as.matrix(expm(G*tsample))


## -----------------------------------------------------------------------------
contour(xc,xc,P/rep(dx,rep(nc,nc)),xlab="From",ylab="To")
abline(0,1,lty="dashed")


## -----------------------------------------------------------------------------
mu <- P %*% xc 
dmudt <- (mu - xc) / tsample
V  <- P %*% (xc^2) - mu^2 
dVdt <- V/tsample

par(mfrow=c(2,1))
plot(xc,dmudt,type="l",lwd=3,xlab="x",ylab="Drift")
lines(xc,f(xc))

plot(xc,dVdt,type="l",lwd=3,xlab="x",ylab="Growth rate in variance")
lines(xc,g(xc)^2)


## -----------------------------------------------------------------------------
obs <- read.table("hmm-obs.txt",header=TRUE)

plot(Y ~t, data=obs)


## -----------------------------------------------------------------------------
vsample <- 0.5
dl <- function(x,y) dpois(y,lambda=vsample*x)
maxy <- max(obs$Y)
ys <- 0:maxy
dltab = outer(xc,ys,dl);

matplot(xc,dltab,type="l")
legend("topright",lty=1+ys,legend=ys,col=1+ys)


## -----------------------------------------------------------------------------
ltab <- dltab[,obs$Y+1]
image(xc,1:length(obs$Y),ltab)  


## -----------------------------------------------------------------------------
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

    return(list(c=const,phi=as,matrix(phi),psi=as.matrix(psi),loglik=sum(log(const))))
}

## Run the filter with default values for parameters    
est <- hmmfilter(G)


## -----------------------------------------------------------------------------
## Show posterior probabilities, and posterior mean 
image(obs$t,xc,est$psi,zlim=c(0,0.2))
lines(obs$t,apply(est$psi,1,function(p)sum(p*xc)))

## Add true state 
sim <- read.table("hmm-states.txt",header=TRUE)
lines(sim$t,sim$X,lwd=2)


## -----------------------------------------------------------------------------
## Likelihood estimation of xi
xis <- seq(0,2,0.05)   # Try these values of xi

loglik <- function(xi)    # Likelihood function of one value of xi 
{
    f <- function(x) lambda*(xi-x);
    D <- function(x) 0.5*gamma^2*x;
    u <- function(x) f(x) - 0.5*gamma^2;


    G = fvade(u,D,xii,'r');

    return(hmmfilter(G)$loglik)
}

## Tabulate likelihood function
ls <- sapply(xis,loglik)

plot(xis,ls,type="b")
grid()
abline(v=xi) ## Include the true value used in the simulation 
abline(h=max(ls)-0.5*qchisq(0.95,df=1)) ## Include confidence interval

