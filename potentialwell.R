# Figure to draw potential well

# The potential
U <- function(x) x^2 + 1/(1-x^2)

# Range of diffusivities (corr. to inverse temperatures)

Ds <- c(0.1,1,10,100)

xmin <- -0.9
xmax <- 0.9


par(mfrow=c(2,1))

plot(U,from=xmin,to=xmax)

Xgrid <- seq(xmin,xmax,length=1000)

phi <- exp(-U(Xgrid)/Ds[1])
phi <- phi / sum(phi) / (xmax-xmin) * length(phi)

plot(Xgrid,phi,type="l",xlab="x",ylab=expression(phi(x)))

print(max(phi))

for(D in Ds[-1])
{
  phi <- exp(-U(Xgrid)/D) / length(Xgrid) / (xmax-xmin)

  phi <- phi / sum(phi) / (xmax-xmin) * length(phi)

  lines(Xgrid,phi)

  print(max(phi))

}

