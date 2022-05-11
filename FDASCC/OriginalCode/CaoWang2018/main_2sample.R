##############################
# confidence band
##############################

rm(list=ls(all=TRUE))

library(MASS)
library(splines)

### Include source files.
source("corr.R")
source("xknots2.R")
source("Ghat.R")
source("Ghat_Nonpara2.R")
source("generateQ_Non_2sample2.R")
source("generateQ2sample2.R")
source("generateCovM.R")
source("cov_rate_2sample.R")
source("cov_rate2sample.R")
source("var_zeta.R")
source("rho_est.R")
source("data_generator2sample.R")
########################
## initial setting
########################
nsample=c(100,100)
n_num=matrix(nsample,2,1,byrow=T)
N_num=n_num						  ## rep(40,length(n_num))
J=3  ### repeated measurements number
select=F

sigma=0.1
Kappa=kappa=4
percent=0.95
mpr=0.9
structure=c("TOEP") #"IND", "EX", "AR1",

structure.true=c("TOEP") #"IND", "EX", "AR1",
if(structure.true=="IND") {rho.true=matrix(rep(1,(J-1)) , 1,(J-1), byrow=T) }
if(structure.true=="TOEP") {rho.true=matrix(c(0.1^seq(1,(J-1))+0.05), 1,(J-1), byrow=T) }
if(structure.true=="AR1") {rho.true=matrix(c(0.2),1,1)}

delta.true=seq(0,0.24,0.03)
set.seed(1000)
Order1=Order2=3
alpha=c(0.05)
nsims=500							 ## no. of replications for coverage frequence

boot=1000	### generate supermum

Ntest=60							## test points
covrate=array(0,c(length(delta.true),nsims,2))
cov.rate=matrix(0,length(delta.true),2)

##############################################
####			start simulation		####
##############################################

for(kkk in 1:length(n_num[1,])){	
	n1=n_num[1,kkk]	### sample size
  n2=n_num[2,kkk]
	N=Ntest=N_num[1,kkk]	#### Xij, j=1,...,N
		
	#phiavgY=matrix(0,N,kappa)
	#phiavg=matrix(0,Ntest,kappa)	 

	#############################################
	#####			knots matrix		###
	#############################################
	 
	knots=knots2(N,Ntest,order1,order2,n1,J)
	XB1=knots$XB1
	beta=knots$beta
	Beta=knots$Beta
	xknots=knots$xknots
	xtestknots=knots$xtestknots
  XBtest1=knots$XBtest1

	#############################################
	###			generate data		###
	#############################################
 for(rr in 1:length(delta.true))
    {  delta=delta.true[rr]
      for(kk in 1:nsims){
			data=data_generator(n1,J,N,Ntest,Kappa,sigma,rho.true,structure,delta,mpr)
			r1=data$select
			nm=sum(r1)
			YY1=data$YY
			Y1=matrix(apply(YY1,2,sum)/nm,N,1)

	 		#############################################
	 		###		start estimation mhat		###
	 		#############################################
			beta.est=beta%*%Y1
			mhat.data1=XB1%*%beta.est

			#############################################
			###	start estimation Ghat (jj1,jj2)	###
			#############################################
			Gest=GhatNP(YY1,mhat.data1,Beta,xknots,xtestknots,Ntest,r1)
			lamdaK.data1=Gest$lamdaKY
			phiavg.data1=Gest$phiavgY
			covM1=Gest$covM 
			GMean1=Gest$GMean      
			kappa.select1=Gest$kappa.select
			num1=Gest$num
			
      data=data_generator(n1,J,N,Ntest,Kappa,sigma,rho.true,structure,0,mpr)
      r2=data$select
      nm=sum(r2)
      YY2=data$YY
      Y2=matrix(apply(YY2,2,sum)/nm,N,1)

	 		#############################################
	 		###		start estimation mhat		###
	 		#############################################
			beta.est=beta%*%Y2
			mhat.data2=XB1%*%beta.est
		        
			#############################################
			###	start estimation Ghat (jj1,jj2)	###
			#############################################
			Gest=GhatNP(YY2,mhat.data2,Beta,xknots,xtestknots,Ntest,r2)
			lamdaK.data2=Gest$lamdaKY
			phiavg.data2=Gest$phiavgY
			GMean2=Gest$GMean                 
			covM2=Gest$covM  
			kappa.select2=Gest$kappa.select
			num2=Gest$num
			
			Q.est=generateQNon2sample(J,covM1,phiavg.data1,lamdaK.data1,GMean1, covM2,phiavg.data2,lamdaK.data2, GMean2,alpha,kappa.select1, num1,kappa.select2, num2)
			Q=Q.est$Q
			covrate[rr,kk,1]=cov_rate(mhat.data1,mhat.data2,GMean1,GMean2,Q,alpha,n1,n2,N)

					ST="IND"                           
					rho.est1=rho_est(YY1,ST,r1,Y1)  
					Gest=Ghat(YY1,mhat.data1,Beta,xknots,xtestknots,Ntest,r1)
					lamdaK.data1=Gest$lamdaKY
					phiavg.data1=Gest$phiavgY
					GMean1=Gest$GMean                 
			 
					kappa.select1=Gest$kappa.num[1]
					num1=Gest$kappa.num[2]
					
					var.zeta1=var_zeta(phiavg.data1,lamdaK.data1,rho.est1,J,ST)
					Gest=Ghat(YY2,mhat.data2,Beta,xknots,xtestknots,Ntest,r2)
					rho.est2=rho_est(YY2,ST,r2,Y2)  
					lamdaK.data2=Gest$lamdaKY
					phiavg.data2=Gest$phiavgY
					var.zeta2=var_zeta(phiavg.data2,lamdaK.data2,rho.est2,J,ST)
					GMean2=Gest$GMean                 
		
					kappa.select2=Gest$kappa.num[1]
					num2=Gest$kappa.num[2]
					
          covM1=generateCovM(J,rho.est1,ST,phiavg.data1,lamdaK.data1,var.zeta1)$covM
          covM2=generateCovM(J,rho.est2,ST,phiavg.data2,lamdaK.data2,var.zeta2)$covM  
	        Q.est=generateQ(covM1, covM2,var.zeta1,var.zeta2,lamdaK.data1,lamdaK.data2,phiavg.data1,phiavg.data2, n1,n2,alpha,kappa.select1, num1,kappa.select2, num2)
					Q=Q.est$Q

					covrate[rr,kk,2]=cov_rate(mhat.data1,mhat.data2,var.zeta1,var.zeta2,Q,alpha,n1,n2,N)
				  
			}
               cov.rate[rr,]=colMeans(covrate[rr,,]) 
               cat("&", (cov.rate[rr,1]))
               cat("&", (cov.rate[rr,2]))
		}
      }
                
                cat("\n")
               for(rr in 1:length(delta.true))  { cat("&", (cov.rate[rr,1])) } 
                  cat("\n")
               for(rr in 1:length(delta.true))  {    cat("&", (cov.rate[rr,2])) } 
 
        write.csv( cov.rate, paste("cov100",structure,".csv", sep="",collapse=NULL))
                  











