
tvec <- seq(10,10000)^4

dt <- diff(tvec)

dB <- rnorm(length(dt),sd=sqrt(dt))

B <- c(0,cumsum(dB)) + rnorm(1,sd=sqrt(tvec[1]))

plot(tvec,B/2/sqrt(tvec*log(log(tvec))),type="l")
