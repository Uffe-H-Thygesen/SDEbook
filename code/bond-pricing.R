require(SDEtools)

## Compute the present-day value of a zero-coupon bond using the Feynman-Kac formula
## assuming risk neutrality, no default, CIR process for the spot interest rate

## The CIR model 
f <- function(x) lambda*(xi-x)
g <- function(x) gamma*sqrt(x)

## Parameters
xi <- 0.025    ## Mean interest rate
lambda <- 0.1 ## 1/(time scale) for fluctuations in the rate
gamma <- 0.00005 ## Noise parameters control fluctuations in rate

## Grid 
xv <- seq(0,0.6,length=101)^3
xc <- 0.5*(head(xv,-1)+tail(xv,-1))


## Times to maturity
tvec <- seq(0,10,1/12)

## Generator
G <- fvade(f,g,xv,'r')

## Plot stationary distribution of the interest rate
require(MASS)
pi <- Null(G)
pi <- pi /diff(xv) / sum(pi)  
plot(xc,pi)


## Sub-generator
Gd <- G - Diagonal(nrow(G),x=xc)

## One time step
dt <- diff(tvec[1:2])
Gdt <- expm(Gd*dt)

## Solve the Feynman-Kac equation 
U <- array(NA,c(nrow(G),length(tvec)))
U[,length(tvec)] <- 1
for(i in 2:length(tvec))
    U[,length(tvec)-i+1] <- as.numeric(Gdt %*% U[,length(tvec)-i+2])

## Plot only in this window
I <- xc < 0.1

pdf(file="bond-pricing.pdf",width=5,height=4)
contour(-rev(tvec),xc[I],t(U[I,]),xlab="Time relative to maturity",ylab="Current spot rate")
dev.off()

## Analytical solution (Sinkala, 2008); compare also (Hull, 2013)
if(FALSE)
{
    gam <- sqrt(lambda^2 + 2 * gamma^2)
    tv <- rev(tvec)
    den <- (gam+lambda)*(exp(gam*tv)-1)+2*gam
    a <- ((2*gam*exp((gam+ lambda)*tv/2)) / den )^(2*lambda*xi/gamma^2)
    b <- 2*(exp(gam*tv)-1)/ den 
    Ua <- t(a*t(exp(-outer(xc,b))))
    contour(-tv,xc[I],t(Ua[I,]),add=TRUE,col="red")
}

