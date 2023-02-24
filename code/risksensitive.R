### Script to illustrate dynamic programming: Risk sensitive investment

### State space: 1..N, the value of a portfolio
### Control: Choose between a safe investment and a risky one.

### Objective: Maximize expected utility


### Number of states
N <-100

### Number of time steps (years)
T <- 500

### Utility function
h <- function(x) (x) # (log(x)+1)

### Mean and mean square of gains in risky and safe strategies

EVtab <-list(c(0.2,0.3),c(0.6,0.8))

Moments2probs <- function(EX,EX2=NULL) 
    {
        if(is.null(EX2))
            {
                EX2 <- EX[2]
                EX <- EX[1]
            }
        
        solve(matrix(c(1,-1,1,1,0,0,1,1,1),c(3,3)),c(1,EX,EX2))
    }

probs2gen <- function(pvec)
    {
        P <- array(0,c(N,N))
        for(i in 1:(N-1))
            {
                P[i,i+1] <- pvec[3]
                P[i+1,i] <- pvec[1]
            }
        P[1,2] <- 0
        P[N,N-1] <- 0
        P <- P + diag(1-apply(P,1,sum))
                
    }

PU <- lapply(EVtab,function(asdf)probs2gen(Moments2probs(asdf)))

print(lapply(PU,max))
print(lapply(PU,min))

### Backward

U <- V <- array(NA,c(N,T))
V[,T] <- h(1:N)

for(t in (T-1):1)
    {
        VU <- do.call("cbind",lapply(PU,function(P) P %*% V[,t+1]))
        U[,t] <- apply(VU,1,which.max)
        V[,t] <- apply(VU,1,max)
    }

require(fields)

par(mfrow=c(2,1))
image.plot(1:T,1:N,t(U))
image.plot(1:T,1:N,t(V))

