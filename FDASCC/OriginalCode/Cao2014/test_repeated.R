
rm(list=ls(all=TRUE))
source("data_generator.R")
source("RepeatedFDAV.R")
source("corr.R")


######################################################################
#    Start generating/reading YY (observed functional data Y(t) )     #
######################################################################
n=50
N=30
J=10
Kappa=2
sigma=0.1
rho.true=0.1
structure.True= "EX"   ## selected structure pool ("IND", "AR1", "EX","TOEP")
X=seq(1/N,1,1/N)
data=data_generator(n,J,N,N,Kappa,sigma,rho.true,structure.True)
YY=data$YY  
##YY is the observed data and the format of it is a J*N*n array, J is the repeated observed curved for each subject. 
#Total number of subject is n and each subject is repeatedly observed J times at N locations (1/N,2/N,..,N/N). So each subject is obersved with a J*N matrix
######################################################################
#   results    #
######################################################################
RESCB(YY,0.05)
