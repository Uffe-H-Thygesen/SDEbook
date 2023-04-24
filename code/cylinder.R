### Water flowing and particles diffusing past a cylinder
### 
### Generates a figure containing:
### . Streamlines past the cylinder (potential flow)
### . A Monte Carlo trajectory of a diffusion particle

### Parameters

Xmax <- 5           ## Plotting region
Ymax <- 3           

### Initialize plot with cylinder / circle
pdf("cylinder.pdf",width=7,height=4)
par(mar=c(4,3,3,1)+0.1)
plot(c(-Xmax,Xmax),c(-Ymax,Ymax),type="n",asp=1,xlab="",ylab="")
phi <- seq(0,2*pi,length=101)
lines(cos(phi),sin(phi),type="l",lwd=5)

### Flow field in polar coordinates
drdt <- function(r,theta) -(1-1/r^2)*cos(theta)
dthetadt <- function(r,theta) (1/r+1/r^3)*sin(theta)

flow <- function(r,theta) c(drdt(r,theta),dthetadt(r,theta))

require(deSolve)

### Plot streamlines by solving ODE's for particle motion
### Note: These could also be derived from the stream function, but such a plot
### would require more computations to look nice and smooth
x0 <- -2*Xmax
for(y0 in seq(0,Ymax,0.5))
    {
        ## Solve for the trajectory in polar coordinates
        sol <- lsoda(c(sqrt(x0^2+y0^2),atan(y0/x0)),seq(0,50,0.01),
                   function(t,y,p) list(flow(y[1],y[2]),numeric(0)),numeric(0))

        ## Convert to Cartesian and plot using symmetries
        lines( sol[,2]*cos(sol[,3]), sol[,2]*sin(sol[,3]))
        lines(-sol[,2]*cos(sol[,3]), sol[,2]*sin(sol[,3]))
        lines( sol[,2]*cos(sol[,3]),-sol[,2]*sin(sol[,3]))
        lines(-sol[,2]*cos(sol[,3]),-sol[,2]*sin(sol[,3]))
    }

### Function for simulating path of particle
require(SDEtools)

## Drift
f <- function(x) flow(x[1],x[2])

## Noise term. Note that noise on the angle is stronger, the smaller the radius is
g <- function(x) sigma*diag(c(1,1/x[1]))

## Projection at each time step - reflect if the particle hits the circle
## with radius 1
p <- function(x)pmax(x,c(1,-Inf))

### Simulate particle path 
x0 <- c(r=9,theta=0.01)
tv <- seq(0,100,0.01)
sigma <- 0.1
set.seed(12345)
path <- as.data.frame(SDEtools::euler(f,g,tv,x0,p=p)$X)

## Add to plot
lines(path$r * cos(path$theta), path$r * sin(path$theta),lwd=3)

dev.off()
