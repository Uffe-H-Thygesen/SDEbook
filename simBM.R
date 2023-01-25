test.rBM <- function(N=1e4,t=c(0,0.5,1.5,2))
{
    B <- sapply(1:N,function(i)rBM(t))
    print("Theoretical covariance:")
    print(sapply(t,function(s)pmin(s,t)))
    print("Empirical covariance:")
    print(cov(t(B)))

    %% QQ-plot
    qqnorm(B[length(t)/sqrt(tail(t,1)),])
}


