##############################
# confidence band
#  11/11/2015
# Use pure Nonparametric method estimate G and SCB
##############################

rm(list=ls(all=TRUE))

library(MASS)
library(splines)

### Include source files.
source("corr.R")
source("xknots.R")
source("data_generator.R")
source("Ghat_Nonpara.R")
source("Ghat.R")
source("generateQ_Non.R")
source("generateQ.R")
source("cov_rate.R")
source("var_zeta.R")
source("rho_est.R")
########################
## initial setting
########################
n_num=c(30,50,100)
N_num=n_num  
J=10  			  ### rep(40,length(n_num))
select=F

sigma=0.1
Kappa=kappa=2

structure.true=c("TOEP") #"IND", "EX", "AR1",
if(structure.true=="IND") {rho.true=matrix(rep(1,(J-1)) , 1,(J-1), byrow=T) }
if(structure.true=="TOEP") {rho.true=matrix(c(0.1^seq(1,(J-1))+0.05), 1,(J-1), byrow=T) }
if(structure.true=="BAND")  {rho.true=matrix(c(0.2),1,1)}
if(structure.true=="AR1") {rho.true=matrix(c(0.2),1,1)}
#if(structure.true=="EX") {rho.true=matrix(c(0.05),1,1)}
set.seed(1000)

Order1=Order2=3
alpha=c(0.05,0.01)
nsims=500							 ## no. of replications for coverage frequence
Ntest=100							## test points
mpr=0.9 #seq(1,0.7,by=-0.3/J)        ## 1-mpr: missing probability


struct=c("IND",structure.true)
covrate=array(0,c(length(struct),nsims,length(alpha)))
cov.rate=matrix(0,length(struct),length(alpha))
covrateNP=array(0,c(nsims,length(alpha)))
##############################################
####			start simulation		####
##############################################

for(kkk in 1:length(n_num)){	
  n=n_num[kkk]	### sample size
  N=N_num[kkk]	#### Xij, j=1,...,N  
  phiavgY=matrix(0,N,kappa)
  phiavg=matrix(0,Ntest,kappa)	 
  
  #############################################
  #####			knots matrix		###
  #############################################
  
  knots=knots2(N,Ntest,order1,order2,n,J)
  XB1=knots$XB1
  beta=knots$beta
  Beta=knots$Beta
  xknots=knots$xknots
  xtestknots=knots$xtestknots
  XBtest1=knots$XBtest1
  
  #############################################
  ###			generate data		###
  #############################################  
  
  
  for(kk in 1:nsims){
    data=data_generator(n,J,N,Ntest,Kappa,sigma,rho.true,structure.true,mpr)
    YY=data$YY
    r=data$select
    nm=sum(r)
    YY=data$YY
    Mtest=data$Mtest
    Y=matrix(apply(YY,2,sum)/nm,N,1)
    
    #############################################
    ###		start estimation mhat		###
    #############################################
    beta.est=beta%*%Y
    mhat.data=XB1%*%beta.est
    mhat.test=XBtest1%*%beta.est
    #############################################
    ###	start estimation Ghat (jj1,jj2)	###
    #############################################

    #print(mean(GMean))
    for(str in 1:length(struct)){
      Gest=Ghat(YY,mhat.data,Beta,xknots,xtestknots,Ntest,r)
      phiavgY.data=Gest$phiavgY
      lamdaK.data=Gest$lamdaK
      phiavg.data=Gest$phiavg
      ST=struct[str]                             
      rho.est=rho_est(YY,ST,r,Y)                          
      var.zeta=var_zeta(phiavg.data,lamdaK.data,rho.est,J,ST)
      #print(mean(var.zeta))
      Q.est=generateQ(J,rho.est,ST,phiavg.data,lamdaK.data,var.zeta,alpha)
      Q=Q.est$Q
      covM=Q.est$covM
      covrate[str,kk,]=cov_rate(mhat.test,Mtest,var.zeta,Q,alpha,n)
    }
    
    for(str in 1:length(struct))  {cov.rate[str,]=colMeans(covrate[str,,]) }
  }
  
  cat("n=", n,"\n")
  cat("structure=", struct[str],"\n")
  for(str in 1:length(struct)) { cat("&", cov.rate[str,1]) } 
  cat("\n")
  for(str in 1:length(struct)) {cat("&", cov.rate[str,2])}
  cat("\n")

}







