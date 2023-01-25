# For a general scalar diffusion on a bounded domain, this script computes the
# probability of exiting to the right, and the expected time to exit.
#
# This is done by formulating the elliptic boundary value problem
# as an initial value problem, applying the standard ODE solver,
# and manipulate the solution to satisfy boundary conditions.

# Define diffusion
f <- function(x) 1
g <- function(x) sqrt(2)

# Interval (a,b)
a <- 10
b <- 20

require(deSolve)

# Write up the second order system 0.5*g^2*s'' + g*s' = 0
# in first order vector form (s,s')
Pdynsys <- function(x,sr,par) list(c(sr[2],-2*sr[2]*f(x)/g(x)/g(x)),numeric(0))

# Discretization
nx <- 101

# Solve for scale function 
s <- lsoda(c(0,1),seq(a,b,length=nx),Pdynsys)

# Extract solution and normalize 
Xvec <- s[,1]
Svec <- s[,2]/s[nx,2]
plot(Xvec,Svec,type="l",xlab="x",ylab=expression({P^x}(X[tau]==b)))

# Now for the expected time to exit
# The system 0.5*g^2*h'' + f*h' + 1 = 0 in first order vector form (h,h')
Edynsys <- function(x,hhp,par) list(c(hhp[2],-2*(hhp[2]*f(x)+1)/g(x)/g(x)),numeric(0))

# Solve for particular solution on *same* grid
hsol <- lsoda(c(0,1),Xvec,Edynsys)

# Extract solution.
# This particular solution satisfies the ODE and the left boundary condition,
# but not (necessarily) the right.
# So:
# Subtract suitable multiple of s(x) to satisfy also right boundary condition
# h(Xmax)=0

Evec <- hsol[,2]-Svec*hsol[nx,2]
plot(Xvec,Evec,type="l",xlab="x",ylab=expression({E^x}*tau))


