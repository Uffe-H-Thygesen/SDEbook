###
### Figure for Gauss curves, illustration for "transport" chapter.
###
###

D <- 1
tvec <- c(1,5,10)

xmin <- -5
xmax <- 5

plot(function(x)dnorm(x,sd=sqrt(2*D*tvec[1])),from=xmin,to=xmax)

for(i in 2:length(tvec))
  plot(function(x)dnorm(x,sd=sqrt(2*D*tvec[i])),from=xmin,to=xmax,add=T)
