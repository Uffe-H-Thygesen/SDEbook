### Monte Carlo analysis of maximum of BM

######################################
### Simulation parameters
######################################

### No Sample paths
M <- 2000

### Duration (no. time steps)
N <- 1000

### Sample time step
dt <- 1

######################################
### Simulate tracks
######################################

### Simulate increments
dB <- array(sqrt(dt)*rnorm(M*N),c(N,M))

### Compute BM track
B <- apply(dB,2,function(db)c(0,cumsum(db)))

### Statistics of each track
S <- apply(B,2,function(B)max(B))
maxAbsB <- apply(B,2,function(B)max(abs(B)))

######################################
### Plot empirical and theoretical
### survivial functions
######################################

### Emprical survival function of S

graphics.off()
plot(sort(S),1-(1:M)/M,type="s",ylim=c(0,1))
### Add theoretical survival function

gridS <- seq(0,max(S),length=1001)
lines(gridS,2*(1-pnorm(gridS/sqrt(N*dt))))

### Emperical survival function max(abs(B))
X11()

plot(sort(maxAbsB),1-(1:M)/M,type="s",ylim=c(0,1),xlim=c(0,max(maxAbsB)))
### Add theoretical survival function
 
gridS <- seq(0,max(S),length=1001)
lines(gridS,4*(1-pnorm(gridS/sqrt(N*dt))),col="blue")

### Yet another approximation. The last  term is abound on the prop of hitting both x anx -x, using reflerction and symmetry, and assuming that we hit x at time tau at t=0 (!) when computing the probability of the transition +x -> -x

PabsXgx <- 4*(1-pnorm(gridS/sqrt(N*dt)))-2*(1-pnorm(gridS/sqrt(N*dt)))*2*(1-pnorm(2*gridS/sqrt(N*dt)))
lines(gridS,PabsXgx,col="red")

