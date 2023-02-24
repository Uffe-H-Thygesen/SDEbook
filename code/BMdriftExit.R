# Figure of probability of exit to the right

pdf(file="BMdriftExit.pdf",width=8,height=4)

par(mfrow=c(1,2),mar=c(5,5,2,2))

Pes <- rev(c(1,10,100))

P <- function(x) (1 - exp(-x*Pe))/(1-exp(-Pe))

Pe <- Pes[1]

plot(P,from=0,to=1,main="",xlab="x/l",ylab=expression({P^x}(X[tau]==l)))

for(i in 2:length(Pes))
{
    Pe <- Pes[i]
    plot(P,from=0,to=1,add=TRUE,lty=i)
}

legend("bottomright",legend=paste("Pe=",Pes),lty=1:length(Pes))

## dev.copy2pdf(file="BMdriftExitP.pdf")

E <- function(x)
     ( (1-x) -  ( exp(- x*Pe ) - exp(-Pe) ) / (1 - exp(-Pe)) ) 

Pe <- Pes[1]

plot(E,from=0,to=1,main="",xlab="x/l",ylab=expression(E^x * tau~" ["*l/u*"]" ))

for(i in 2:length(Pes))
{
    Pe <- Pes[i]
    plot(E,from=0,to=1,add=TRUE,lty=i)
}

## dev.copy2pdf(file="BMdriftExitE.pdf")

dev.off()
