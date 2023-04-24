### Illustration of Gaussian bell attenuated by diffusion


pdf(file="figGauss.pdf",width=7,height=4)

D <- 1                 ## Diffusivity 
tvec <- c(1,5,10)      ## Vector of time points 

xmin <- -10            ## Plotting window
xmax <- -xmin

ltys <- c("solid","dashed","dotted")

## Plot p.d.f. at the first time point
plot(function(x)dnorm(x,sd=sqrt(2*D*tvec[1])),from=xmin,to=xmax,
     xlab="Position [m]",lty=ltys[1],ylab="Concentration [1/m]")

## Add curves for the remaining time points
for(i in 2:length(tvec))
  plot(function(x)dnorm(x,sd=sqrt(2*D*tvec[i])),lty=ltys[i],
       from=xmin,to=xmax,add=TRUE)


legend(4.5,0.25,legend=paste("Time =",tvec,"s"),lty=ltys)

dev.off()

