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
source("Ghat_Nonpara2.R")
source("Ghat.R")
source("generateQ_Non.R")
source("generateQ.R")
source("cov_rate.R")
source("var_zeta.R")
source("rho_est.R")
########################
## initial setting
########################
n_num=c(30, 60,120)
N_num=n_num 
J=5  			  ### rep(40,length(n_num))
select=F
percent=0.95

sigma=0.1
Kappa=kappa=4

structure.true=c("TOEP") #"IND", "EX", "AR1",
if(structure.true=="IND") {rho.true=matrix(rep(1,(J-1)) , 1,(J-1), byrow=T) }
if(structure.true=="TOEP") {rho.true=matrix(c(0.1^seq(1,(J-1))+0.05), 1,(J-1), byrow=T) }
if(structure.true=="BAND")  {rho.true=matrix(c(0.2),1,1)}
if(structure.true=="AR1") {rho.true=matrix(c(0.2),1,1)}
if(structure.true=="EX") {rho.true=matrix(c(0.05),1,1)}
set.seed(2000)

Order1=Order2=3
alpha=c(0.05,0.01)
nsims=500							 ## no. of replications for coverage frequence
#Ntest=100							## test points
mpr=seq(1,0.8,by=-0.2/J)       ## 1-mpr: missing probability


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
  Ntest=N
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
      if (mean(mpr)<1)
      {YYi=matrix(NA, n,N)
        for ( i in 1:n) {YYi[i,]= colMeans(YY[,,i])*J/sum(r[i,])}
        Y= colMeans(YYi)
        }
     
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
      Gest=GhatNP(YY,mhat.data,Beta,xknots,xtestknots,Ntest,r)
      phiavgY.data=Gest$phiavgY
      lamdaK.data=Gest$lamdaK
      phiavg.data=Gest$phiavg
      GMean=Gest$GMean                 
      covM=Gest$covM  
      indexka=Gest$indexka
      kappa.select=Gest$kappa.select
      num=Gest$num
      Q.est=generateQNon(J,covM,phiavg.data,lamdaK.data,GMean,alpha,indexka, kappa.select, num)
      Q=Q.est$Q
      covrateNP[kk,]=cov_rate(mhat.test,Mtest,GMean,Q,alpha,n)   
      for(str in 1:length(struct)){
        Gest=Ghat(YY,mhat.data,Beta,xknots,xtestknots,Ntest,r)
        phiavgY.data=Gest$phiavgY
        lamdaK.data=Gest$lamdaK
        phiavg.data=Gest$phiavg
        kappa.select=(Gest$kappa.num)[1]
        num=(Gest$kappa.num)[2]
        ST=struct[str]                             
        rho.est=rho_est(YY,ST,r,Y)                          
        var.zeta=var_zeta(phiavg.data,lamdaK.data,rho.est,J,ST)
        Q.est=generateQ(J,rho.est,ST,phiavg.data,lamdaK.data,var.zeta,alpha, kappa.select, num)
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
  cov.rateNP=colMeans(covrateNP) 
  cat("Nonparameteric", "\n")
  cat("&", (cov.rateNP[1]))  
  cat("\n")
  cat("&", (cov.rateNP[2]))
  cat("\n")
}







