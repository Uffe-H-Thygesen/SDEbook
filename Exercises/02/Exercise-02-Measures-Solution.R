## -----------------------------------------------------------------------------
N <- 1e4

X <- runif(N)
Y <- runif(N)
Z <- X+Y

dat <- data.frame(X=X,Y=Y,Z=Z)

Xbreaks <- seq(0,1,0.1)
Xc <- Xbreaks[-1] - 0.5*diff(Xbreaks)

EZgX <- tapply(Z,cut(X,breaks=Xbreaks),mean)

plot(Xc,EZgX)
plot(function(x)x+0.5,add=TRUE,from=0,to=1)


## -----------------------------------------------------------------------------
Zbreaks <- seq(0,2,0.1)
Zc <- Zbreaks[-1] - 0.5*diff(Zbreaks)

cZ <- cut(Z,breaks=Zbreaks)

EXgZ <- tapply(X,cZ,mean)
plot(Zc,EXgZ)
plot(function(z)0.5*z,add=TRUE,from=0,to=2)


## -----------------------------------------------------------------------------
## Based on analytical conditioning:
mean(X)
mean(Z/2)


## -----------------------------------------------------------------------------
## Based on binned Z's:
mean(X)
pZ <- table(cZ) / length(Z)
sum(EXgZ*pZ)


## -----------------------------------------------------------------------------
## Conditional variance
VXgZ <- tapply(X,cZ,var)

plot(Zc,VXgZ)
VXgzanalytical <- function(z) pmin(z,2-z)^2/12
curve(VXgzanalytical,add=TRUE)


## -----------------------------------------------------------------------------
## Variance decomposition, analytically:
print(var(X))
print(mean(VXgzanalytical(Z)) + var(Z/2))

## Variance decomposition, purely empirical:
VXgZ[is.na(pZ)] <- 0
EXgZ[is.na(pZ)] <- 0
pZ[is.na(pZ)] <- 0

print(sum(pZ*VXgZ) + sum(pZ*(EXgZ-sum(pZ*EXgZ))^2))

