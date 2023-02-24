require(SDEtools)

T <- 100
dt <- 2^{-10}

tvec <- seq(0,T,dt)

xi <- 1
lambda <- 1
gamma <- 0.1

x0 <- 0.1

f <- function(x) lambda*(xi-x)
g <- function(x) gamma*sqrt(x)

B <- rBM(tvec)

sim <- euler(f,g,x0,tvec,B)

plot(sim)
