### LQR Tracking a quadratic reference

rm(list=ls())
graphics.off()


### System dynamics
a <- -1
b <- 1

### Reference dynamics
aR <- 0

A <- matrix(c(0,0,0,aR),nrow=2,ncol=2)
F <- matrix(c(1,0),nrow=2,ncol=1)
G <- matrix(c(0,1),nrow=2,ncol=2)

Q <- matrix(c(1,-1,-1,1),nrow=2,ncol=2)
R <- 1
P <- NULL

nx <- 3
nu <- 2
nb <- 3


A <- matrix(rnorm(nx*nx),nrow=nx,ncol=nx)
F <- matrix(rnorm(nx*nu),nrow=nx,ncol=nu)
G <- matrix(rnorm(nx*nb),nrow=nx,ncol=nb)

Q <- matrix(rnorm(nx*nx),nrow=nx,ncol=nx)
Q <- Q %*% t(Q)

R <- matrix(rnorm(nu*nu),nrow=nu,ncol=nu)
R <- R %*% t(R)

P <- matrix(rnorm(nx*nx),nrow=nx,ncol=nx)
P <- P %*% t(P)

Riccati <- function(time,Ss,par)
    {
        S <- matrix(Ss[1:(length(Ss)-1)],nrow=nrow(par$A))
        SA <- S %*% par$A
        dSdt <- SA + t(SA) - S %*% par$FRiF %*% S + par$Q
        dsdt <- 0.5*sum(diag(par$Gt %*% S %*% par$G))

        return(list(c(as.numeric(dSdt),dsdt)))
    }


lqr <- function(times,A,F,G,Q,R,P=NULL)
{
        RiF <- solve(R,t(F))

        if(is.null(P)) P <- array(0,rep(nrow(A),2))
        
        par <- list(A = A,G = G,Gt = t(G),FRiF = F %*% RiF,Q=Q)
        
        require(deSolve)

        sol <- ode(y = c(as.numeric(P),0),times=times,func=Riccati,parms=par)

        SvecToGain <- function(Ss) as.numeric(-RiF %*% matrix(Ss[2:(length(Ss)-1)],nrow=nrow(par$A)))

        s <- sol[,ncol(sol)]
        S <- array(sol[,2:(ncol(sol)-1)],c(length(times),nrow(A),nrow(A)))
        L <- apply(sol,1,SvecToGain)

        L <- array(t(L),c(length(times),ncol(F),nrow(A)))

        return(list(times=times,A=A,F=F,G=G,Q=Q,R=R,P=P,S=S,s=s,L=L))

    }

times <- seq(0,30,0.01)/10
sol <- lqr(times,A,F,G,Q,R,P)
        
# matplot(sol$times, [,1],sol[,1+(1:length(A))],type="l")

matplot(sol$times,array(sol$L,c(length(times),2)),type="l")

plot(sol$times,sol$L[,1,2],type="l")

### Beam example

L <- 10
N <- 10
dx <- L/N
EI <- 1
rhoA <- 1

e1 <- 0.1
nu <- 0.1
    
dx2 <- (L/N)^2

mydiag <- function(x,k)
    {
        if(k==0) return(diag(x))
        if(k>0) return(cbind(rbind(array(0,c(k,length(x))),diag(x)),array(0,c(length(x)+k,k))))
        if(k<0) return(t(mydiag(x,-k)))
    }

y2M <- EI/dx2*rbind(c(2,rep(0,N-1)),mydiag(rep(1,N-1),-1) + mydiag(c(rep(1,N-2),0),1) - diag(c(2,rep(2,N-2),0)))

M2a <- 1/rhoA/dx2*cbind(c(1,rep(0,N-1)),diag(rep(-2,N)) + mydiag(rep(1,N-1),-1) + mydiag(c(rep(1,N-2,),2),1))

K <- M2a %*% y2M
evs <- eigen(K)
matplot(evs$vectors[,N:(N-4)],type="l")

O <- array(0,c(N,N))
I <- diag(rep(1,N))
A <- rbind(cbind(O,I),
           cbind(-K,-e1*K-nu*diag(rep(1,N))))

Q <- cbind(y2M,array(0,c(N+1,N)))
Q <- t(Q) %*% Q

Q <- diag(rep(1,nrow(A)))

F <- matrix(c(rep(0,N),M2a[,1]),nrow=2*N,ncol=1)
G <- matrix(0,nrow=2*N,ncol=1)
G[2*N,1] <- 1
R <- 1
P <-  NULL

times <- seq(0,100,1)
sol <- lqr(times,A,F,G,Q,R,P)
        
# matplot(sol$times, [,1],sol[,1+(1:length(A))],type="l")

matplot(sol$times,array(sol$L,c(length(times),2*N)),type="l")

X11()
par(mfrow=c(2,1))
plot(sol$L[length(times),1,1:N])
plot(sol$L[length(times),1,N+(1:N)])
