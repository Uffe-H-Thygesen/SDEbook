p <- function(s,t,x,y)
  {
    mu <- exp(-lambda*(t-s))*x
    Sigma <- sigma^2/2/lambda*(1-exp(-2*lambda*(t-s)))

    return(dnorm(y,mean=mu,sd=sqrt(Sigma)))
  }

graphics.off()

Xmin <- -10
Xmax <- 10
Ymin <- -4
Ymax <- 4

lambda <- 1
sigma <- 1
x <- - 3
s <- 0
t <- 0.1
y <- -1

lags <- 0.1*2^seq(0,10)
lags <- lags[-1]

plot(function(y)p(s,lags[1],x,y),
     from=Ymin,to=Ymax,main="Forward transition probabilities",
     n = 1001,
     xlab="Terminal position y",
     ylab="Transition density p(0->t,x->y)"
     )

for(i in 2:length(lags))
  {
    plot(function(y)p(0,lags[i],x,y),from=Ymin,to=Ymax,lty=i,add=TRUE)
  }

dev.copy2pdf(file="ou-kolmo-fwd.pdf")

X11()

plot(function(x)p(0,lags[1],x,y),
     from=Xmin,to=Xmax,main="Backward transition probabilities",
     n = 1001,
     xlab="Initial position x",
     ylab="Transition density p(s->0,x->y)"
     )

for(i in 2:length(lags))
  {
    plot(function(x)p(0,lags[i],x,y),from=Xmin,to=Xmax,lty=i,add=TRUE)
  }

dev.copy2pdf(file="ou-kolmo-bkwd.pdf")
