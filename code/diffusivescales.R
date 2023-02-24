D = c(1e-9,2e-5,2e-11)

T = c(1,60,3600,3600*24)

t(sapply(D,function(DD)sqrt(2*DD*T)))

