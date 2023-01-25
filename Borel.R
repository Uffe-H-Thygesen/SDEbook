### Illustration of Borel's paradox.

graphics.off()
par(cex=1.5)

xv <- seq(-2,2,0.05)

phi <- dnorm(xv)

Phi <- outer(phi,phi)


col <- rgb(rev(0:255),255,rev(0:255),255,max=255)

image(xv,xv,Phi,xaxt="n",yaxt="n",col=col,xlab="",ylab="",asp=1,bty="n")


### Insert contour lines
### X^2+Y^2 is chisq(2) distributed, so \P(S <= s) = dchisq

Squantiles <- c(0.5)

thetas <- seq(0,2*pi,length=101)

for(s in qchisq(Squantiles,df=2))
  lines(sqrt(s)*cos(thetas),sqrt(s)*sin(thetas),col="grey")

lines(xv,numeric(length(xv)),lwd=5)


X <- 1.5
Y <- 0.8

points(X,Y,pch=16,cex=1)

text(X,Y,expression("("*X*","*Y*")"),pos=3)

arrowhead <- function(x,y,angle,size,sharpness=0.5,...)
  {
    xs <- size*c(-1,0,-1)
    ys <- size*sharpness*c(1,0,-1)

    polygon(x+xs*cos(angle)-ys*sin(angle),
            y+xs*sin(angle)+ys*cos(angle),...)
  }

lines(c(0,X),c(0,Y),lwd=2)
arrowhead(0.98*X,0.98*Y,atan2(Y,X),0.1,col="black")
          
text(X*0.5,Y*0.5,"R",pos=3)


R <- sqrt(X^2+Y^2)
Theta <- atan2(Y,X)

thetas <- seq(0,Theta,length=20)

lines(R/2*cos(thetas),R/2*sin(thetas))

arrowhead(R/2*cos(Theta),R/2*sin(Theta),Theta+pi/2,size=0.1,col="black")
text(R/2*cos(Theta/2),R/2*sin(Theta/2),expression(Theta),pos=4)

dev.copy2pdf(file="BorelSketch.pdf")

