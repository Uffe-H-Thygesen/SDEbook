norminf <- function(A)
# Return the infinity norm of matrix A, i.e. the largest row sum
# of abs(A)
{
 max(apply(abs(A),1,sum))
}

### Generate a random transition matrix for a
### continuous-time markov chain
rmarkov <- function(n = 10)
{
 A <- matrix(runif(n*n),n,n)
 diag(A) <- 0
 diag(A) <- -apply(A,1,sum)
 A
}

expm.pade <- function(A, q = 6)
# Returns matrix exponential using a scale-and-square algorithm
{
# Rescale
 ni <- norminf(A)
 e <- ceiling(log2(ni))
 f <- ni/2^e
 s <- max(0,e+1)
 A <- A/2^s

# Pade approximation for exp(A)
 X <- A
 c <- 0.5
 E <- diag(rep(1,nrow(A))) + c*A
 D <- diag(rep(1,nrow(A))) - c*A

 p <- TRUE

 for(k in 2:q)
 {
  c <- c* (q - k + 1) / ( k*(2*q - k+1))
  X <- A %*% X
  cX <- c*X
  E <- E + cX
  if(p)
   D <- D + cX
  else
   D <- D - cX

  p <- !p
 }

 E <- solve(D,E)
 for(k in 1:s) E <- E %*% E

 return(E)
}

# Default choice for matrix exponential: The Pade
# algorithm
expm <- expm.pade

# An eigen-decomposition method for the matrix exponential
# Note: This requires that A is diagonizable (i.e., no Jordan blocks)
expm.eigen <- function(A,return.real = is.real(A))
{
 e <- eigen(A)
 e <- t(solve(t(e$vec),t(e$vectors %*% diag(exp(e$values)))))
 if(return.real) e <- Re(e)
 return(e)
}


# Test comparison between the two methods.
# Note: rmarkov() returns a diagonizable matrix (w.p. 1)
test.expm <- function(m = rmarkov())
{
 e1 <- expm.pade(m)
 e2 <- expm.eigen(m)

 norminf(e1 - e2) 
}

test.expm.pade <- function(m = rmarkov(),nmax = 8)
{
 e <- expm.pade(m,nmax)
 is <- 1:(nmax-1)
 es <- sapply(is,function(i)norminf(e-expm.pade(m,i)))
 plot(is,es,log = "y")
 list(is,es)
}

stat.mc <- function(nw = 72,times = seq(0,20,4))
{
 nt <- length(times)
 A <- function(i) expm.pade(rmarkov(nw))

 TA <- diag(rep(-1,nt*nw))
 for(i in 1:nt) 
  TA[1:nw + (i - 1 )*nw,1:nw + (i %% nt)*nw] <- A(i)
 TA
}
