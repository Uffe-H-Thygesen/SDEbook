## Resonance

d <- 1
w0 <- 1

H <- function(w) { s <- 1i*w ; 1/(s^2 + 2*d*w0*w + w0^2) }
S <- function(w) abs(H(w))^2

plot(S,from=1e-2,to=1e2,log="xy")

## State space implementation, compare with the exercise on the noisy oscillator

## Simulate the system

## Repeat for the other ones ...
