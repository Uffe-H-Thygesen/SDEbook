T <- 1000
d <- numeric(T)
l <- 2
mu <- 1
sigma <- 1
for(i in 2:T) d[i] <- l*d[i-1] + mu + rnorm(1,sd=sigma)
plot(d,log="y")
plot(d^2,log="y")

i <- 1:T
lines(i,((i^2*mu^2+2*i*sqrt(i)*mu*sigma + sigma^2*i)*l^(2*i)))
