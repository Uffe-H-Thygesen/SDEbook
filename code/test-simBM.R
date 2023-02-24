
R version 3.6.1 (2019-07-05) -- "Action of the Toes"
Copyright (C) 2019 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> setwd('/home/uhth/Projects/sdebook/figures/')
> rnorm(sd=c(1,2))
Error in rnorm(sd = c(1, 2)) : argument "n" is missing, with no default
> rBM <- function(times) c(0,cumsum(rnorm(sd = sqrt(dt <- diff(times)),n=length(dt))))
> rBM(0:10)
[1]  0.000000 -0.730934
> rBM(c(1,2))
[1]  0.000000 -1.145846
> rBM <- function(t) c(0,cumsum(rnorm(n=length(t)-1,sd = sqrt(diff(t)))))
> rBM(1:10)
 [1]  0.0000000 -0.2282323  0.3699393  0.1685694  1.8919696  3.9389929
 [7]  3.3654971  2.5477495  2.2832518  2.4009791
> rBM <- function(t) cumsum(rnorm(n=length(t),sd = sqrt(diff(c(0,t)))))
> rBM(1:10)
 [1] -0.6400218 -2.0782423 -2.4986264 -3.3880548 -2.7418932 -3.8490337
 [7] -2.8375454 -3.8141521 -3.8782398 -2.9654447
> rBM(0:10)
 [1]  0.0000000  0.9098639  1.1563019  0.7089319 -0.6328042  1.5851683
 [7]  0.8145799  0.5623929  1.5399605  1.1768345  1.4142137
> rBM <- function(t) cumsum(rnorm(n=length(t),sd = sqrt(diff(c(0,t)))))
> 
> test.rBM <- function(N=1e3,tvec = c(0,0.5,1.5,2))
+ {
+     B <- sapply(1:N,function(i)rBM(tvec)$B)
+ 
+     print("Theoretical covariance:")
+     print(sapply(tvec,function(t)pmin(t,tvec)))
+     print("Empirical covariance:")
+     print(cov(t(B)))
+ }
> test.rBM()
Error in rBM(tvec)$B (from #3) : $ operator is invalid for atomic vectors
> + + + + + + + + + + > rBM <- function(t) cumsum(rnorm(n=length(t),sd = sqrt(diff(c(0,t)))))
> test.rBM <- function(N=1e3,tvec = c(0,0.5,1.5,2))
+ {
+     B <- sapply(1:N,function(i)rBM(tvec))
+ 
+     print("Theoretical covariance:")
+     print(sapply(tvec,function(t)pmin(t,tvec)))
+     print("Empirical covariance:")
+     print(cov(t(B)))
+ }
> test.rBM()
[1] "Theoretical covariance:"
     [,1] [,2] [,3] [,4]
[1,]    0  0.0  0.0  0.0
[2,]    0  0.5  0.5  0.5
[3,]    0  0.5  1.5  1.5
[4,]    0  0.5  1.5  2.0
[1] "Empirical covariance:"
     [,1]      [,2]      [,3]      [,4]
[1,]    0 0.0000000 0.0000000 0.0000000
[2,]    0 0.4723318 0.4485142 0.4238854
[3,]    0 0.4485142 1.3961198 1.3564542
[4,]    0 0.4238854 1.3564542 1.8429861
> + + + + + + + + + + > rBM <- function(t) cumsum(rnorm(n=length(t),sd = sqrt(diff(c(0,t)))))
> test.rBM <- function(N=1e4,tvec = c(0,0.5,1.5,2))
+ {
+     B <- sapply(1:N,function(i)rBM(tvec))
+ 
+     print("Theoretical covariance:")
+     print(sapply(tvec,function(t)pmin(t,tvec)))
+     print("Empirical covariance:")
+     print(cov(t(B)))
+ }
> test.rBM()
[1] "Theoretical covariance:"
     [,1] [,2] [,3] [,4]
[1,]    0  0.0  0.0  0.0
[2,]    0  0.5  0.5  0.5
[3,]    0  0.5  1.5  1.5
[4,]    0  0.5  1.5  2.0
[1] "Empirical covariance:"
     [,1]      [,2]      [,3]      [,4]
[1,]    0 0.0000000 0.0000000 0.0000000
[2,]    0 0.4957874 0.4985028 0.4881005
[3,]    0 0.4985028 1.5334471 1.5297263
[4,]    0 0.4881005 1.5297263 2.0184418
> outer(tvec,tvec,function(x,y)min(x,y))
Error in outer(tvec, tvec, function(x, y) min(x, y)) : 
  object 'tvec' not found
> outer(0:4,0.4,function(x,y)pmin(x,y))
     [,1]
[1,]  0.0
[2,]  0.4
[3,]  0.4
[4,]  0.4
[5,]  0.4
> outer(0:4,0.4,function(x,y)pmin(x,y))  C-c C-c
>   C-c C-c
>   itointegral <- function(f,g) c(0,cumsum(head(f,-1)*diff(g)))

  itointegral <- function(f,g) c(0,cumsum(head(f,-1)*diff(g)))
> 
> 
> itointegral(c(1,2,3),c(1,2,3))
[1] 0 1 3
> test.rBM <- function(N=1e3,t=c(0,0.5,1.5,2))
{
    B <- sapply(t,function(t)rBM(t))
    print(sapply(t,function(s)pmin(s,t)))
    print(cov(B))
}



test.rBM <- function(N=1e3,t=c(0,0.5,1.5,2))
+ {
+     B <- sapply(t,function(t)rBM(t))
+     print(sapply(t,function(s)pmin(s,t)))
+     print(cov(B))
+ }
> 
> 
> 
> 
> test.rBM()
     [,1] [,2] [,3] [,4]
[1,]    0  0.0  0.0  0.0
[2,]    0  0.5  0.5  0.5
[3,]    0  0.5  1.5  1.5
[4,]    0  0.5  1.5  2.0
Error in cov(B) : supply both 'x' and 'y' or a matrix-like 'x'
> + + + + + + + > test.rBM <- function(N=1e3,t=c(0,0.5,1.5,2))
+ {
+     B <- sapply(t,function(t)rBM(t))
+     print("Theoretical covariance:")
+     print(sapply(t,function(s)pmin(s,t)))
+     print("Empirical covariance:")
+     print(cov(t(B)))
+ }
> test.rBM()
[1] "Theoretical covariance:"
     [,1] [,2] [,3] [,4]
[1,]    0  0.0  0.0  0.0
[2,]    0  0.5  0.5  0.5
[3,]    0  0.5  1.5  1.5
[4,]    0  0.5  1.5  2.0
[1] "Empirical covariance:"
     [,1] [,2] [,3] [,4]
[1,]   NA   NA   NA   NA
[2,]   NA   NA   NA   NA
[3,]   NA   NA   NA   NA
[4,]   NA   NA   NA   NA
> + + + + + + + > test.rBM <- function(N=1e3,t=c(0,0.5,1.5,2))
+ {
+     B <- sapply(1:N,function(i)rBM(t))
+     print("Theoretical covariance:")
+     print(sapply(t,function(s)pmin(s,t)))
+     print("Empirical covariance:")
+     print(cov(t(B)))
+ }
> test.rBM()
[1] "Theoretical covariance:"
     [,1] [,2] [,3] [,4]
[1,]    0  0.0  0.0  0.0
[2,]    0  0.5  0.5  0.5
[3,]    0  0.5  1.5  1.5
[4,]    0  0.5  1.5  2.0
[1] "Empirical covariance:"
     [,1]      [,2]      [,3]      [,4]
[1,]    0 0.0000000 0.0000000 0.0000000
[2,]    0 0.5043719 0.4859947 0.5121148
[3,]    0 0.4859947 1.4393329 1.4856308
[4,]    0 0.5121148 1.4856308 2.0029957
> + + + + + + + > test.rBM <- function(N=1e4,t=c(0,0.5,1.5,2))
+ {
+     B <- sapply(1:N,function(i)rBM(t))
+     print("Theoretical covariance:")
+     print(sapply(t,function(s)pmin(s,t)))
+     print("Empirical covariance:")
+     print(cov(t(B)))
+ }
> test.rBM()
[1] "Theoretical covariance:"
     [,1] [,2] [,3] [,4]
[1,]    0  0.0  0.0  0.0
[2,]    0  0.5  0.5  0.5
[3,]    0  0.5  1.5  1.5
[4,]    0  0.5  1.5  2.0
[1] "Empirical covariance:"
     [,1]      [,2]      [,3]      [,4]
[1,]    0 0.0000000 0.0000000 0.0000000
[2,]    0 0.5095855 0.5183415 0.5196074
[3,]    0 0.5183415 1.5248947 1.5127316
[4,]    0 0.5196074 1.5127316 2.0001652
> 