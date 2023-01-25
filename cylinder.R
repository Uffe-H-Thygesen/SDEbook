### Water flowing and particles diffusing past a cylinder
### 
### Generates a figure containing:
### . Streamlines past the cylinder (potential flow)
### . A Monte Carlo trajectory of a diffusion particle

### Parameters

Xmax <- 5
Ymax <- 3

### Initialize plot

graphics.off()
pdf("cylinder.pdf",width=7,height=4)
par(mar=c(4,3,3,1)+0.1)
plot(c(-Xmax,Xmax),c(-Ymax,Ymax),type="n",asp=1,xlab="",ylab="")
phi <- seq(0,2*pi,length=101)
lines(cos(phi),sin(phi),type="l",lwd=5)


### Flow field 
drdt <- function(r,theta) -(1-1/r^2)*cos(theta)
dthetadt <- function(r,theta) (1/r+1/r^3)*sin(theta)

flow <- function(r,theta) c(drdt(r,theta),dthetadt(r,theta))

require(deSolve)

### Plot streamlines
x0 <- -2*Xmax
for(y0 in seq(0,Ymax,0.5))
    {
      sol <- lsoda(c(sqrt(x0^2+y0^2),atan(y0/x0)),seq(0,50,0.01),
                   function(t,y,p) list(flow(y[1],y[2]),numeric(0)),numeric(0))

      lines( sol[,2]*cos(sol[,3]), sol[,2]*sin(sol[,3]))
      lines(-sol[,2]*cos(sol[,3]), sol[,2]*sin(sol[,3]))
      lines( sol[,2]*cos(sol[,3]),-sol[,2]*sin(sol[,3]))
      lines(-sol[,2]*cos(sol[,3]),-sol[,2]*sin(sol[,3]))
    }

### Function for simulating path of particle
euler <- function(sigma=1,r0=8,theta0=0.01,t1=100,dt=0.01)
  {
    nt <- round(t1/dt)
    theta <- r <- numeric(nt+1)
    r[1] <- r0
    theta[1] <- theta0
    
    sdt <- sigma*sqrt(dt)

    for(i in 1:nt)
      {
        r[i+1] <- r[i] + drdt(r[i],theta[i]) * dt + sdt*rnorm(1)

        r[i+1] <- abs(r[i+1]-1)+1
        
        theta[i+1] <- theta[i] + dthetadt(r[i],theta[i]) * dt + sdt/r[i]*rnorm(1)
      }

    return(list(r=r,theta=theta))
  }


### Simulate particle path and add to plot
path <- euler(sigma=0.1)
lines(path$r * cos(path$theta), path$r * sin(path$theta),lwd=3)
print(path$theta[length(path$theta)])

dev.off()
