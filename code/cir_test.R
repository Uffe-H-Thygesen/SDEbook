## Test of the CIR transition probabilities

## OU processes
n <- 10
mu <- 2
sigma <- 1

## CIR model parameters
gamma <- 2*sigma
lambda <- 2*mu
xi <- n*sigma^2/lambda

## Simulation
t <- 10
N <- 10000
x0 <- 2

## Simulation with OU-processes
rcir.ou <- function(t,x0,n,mu,sigma)
{
    ## Initial position of each OU-particle
    y0 <- sqrt(x0/n)

    ## Terminal position of each OU-particle
    yt <- rnorm(n,mean=y0*exp(-mu*t),sd=sigma/sqrt(2*mu)*sqrt((1-exp(-2*mu*t))))

    return(sum(yt^2))
}

## Simulation with Euler-Maruyama
rcir.EM <- function(t,x0,gamma,lambda,xi)
{
    f <- function(x) lambda*(xi-x)
    g <- function(x) gamma*sqrt(abs(x))

    tv  <- seq(0,t,length=101)
    
    sim <- SDEtools::euler(f,g,tv,x0,p=abs)
    return(sim$X[length(tv)])
}

X.ou <- sapply(1:N,function(i)rcir.ou(t,x0,n,mu,sigma))
X.EM <- sapply(1:N,function(i)rcir.EM(t,x0,gamma,lambda,xi))

dcir <- function(xt,t,x0,gamma,lambda,xi)
{
    c <- 2*lambda/gamma^2/(1-exp(-lambda*t))
    n <- 4*lambda*xi/gamma^2
    nu <- 2*c*x0*exp(-lambda*t)

    return(c*exp(-(2*c*xt+nu)/2)*(2*c*xt/nu)^(n/4-1/2)*besselI(sqrt(2*c*nu*xt),n/2-1))
}

dcir.chisq <- function(xt,t,x0,gamma,lambda,xi)
{
    c <- 2*lambda/gamma^2/(1-exp(-lambda*t))
    n <- 4*lambda*xi/gamma^2
    nu <- 2*c*x0*exp(-lambda*t)

    return(dchisq(2*c*xt,df=n,ncp=nu)*2*c)
}

par(mfrow=c(2,1))
xmax <- max(c(X.ou,X.EM))
hist(X.ou,breaks=seq(0,xmax,length=50),freq=FALSE)
plot(function(x)dcir.chisq(x,t,x0,gamma,lambda,xi),from=0,to=xmax,add=TRUE)
plot(function(x)dcir(x,t,x0,gamma,lambda,xi),from=0,to=xmax,add=TRUE,col="red")
hist(X.EM,breaks=seq(0,xmax,length=50),freq=FALSE)
plot(function(x)dcir.chisq(x,t,x0,gamma,lambda,xi),from=0,to=xmax,add=TRUE)
plot(function(x)dcir(x,t,x0,gamma,lambda,xi),from=0,to=xmax,add=TRUE,col="red")

## plot(function(x)dcir(x,t,x0,gamma,lambda,xi),from=0,to=xmax)
## plot(function(x)dcir.chisq(x,t,x0,gamma,lambda,xi),from=0,to=xmax)

xv <- seq(0,xmax,length=100000)
pdf <- dcir.chisq(xv,t,x0,gamma,lambda,xi)
sum(pdf[-1])*(xv[2]-xv[1])
