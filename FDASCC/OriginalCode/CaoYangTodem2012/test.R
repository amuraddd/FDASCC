n=50
N=30
source("SCB.R")
YY=matrix(0,n,N)     ### YY matrix is the observed data and the format of it is n*N
######################################################################
#    Start generating/reading YY (observed functional data Y(t) )     #
######################################################################
X=seq(1/N,1,1/N)
for(i in 1:n)
  { YY[i,]=10*sin(X)+(rnorm(N)) }
SCB(YY,0.05)