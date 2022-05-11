##############################
# confidence band
##############################

rm(list=ls(all=TRUE))

library(MASS)
library(splines)

### Include source files.
source("corr.R")
source("xknots2.R")
source("data_generator.R")
source("Ghat.R")
source("generateQ.R")
source("cov_rate.R")
source("var_zeta.R")
source("rho_est.R")
########################
## initial setting
########################
n_num=c(30,50,100)
N_num=n_num						  ## rep(40,length(n_num))
J=10  ### repeated measurements number
select=T

sigma=0.1
Kappa=kappa=2

structure.true= c("IND", "AR1", "EX","TOEP")

#  structure="AR1"
# structure="EX"
#### structure="BAND"
# structure="TOEP"
# structure="IND"



set.seed(100)

Order1=Order2=3
alpha=c(0.05,0.01)
nsims=500							 ## no. of replications for coverage frequence

boot=1000	### generate supermum
bootnum=500 

Ntest=100							## test points

struct=c( "IND","EX", "AR1","TOEP")  #c("IND", "EX", "AR1","TOEP")
conboot=array(0,c(bootnum,2,length(struct)))
choose_str=matrix(0,length(struct),2)
index=matrix(0,length(n_num),nsims)
covrate=array(0,c(length(struct),nsims,length(alpha)))


#if(select==F){struct<-structure}
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
             
	for(rr in 1:length(structure.true)){
          	      
              
                   structure=structure.true[rr]
                  if(structure=="TOEP") {rho.true=0.5*seq(1,(J-1))/(3*J)}
                  if(structure=="IND") {rho.true=0}
                  if(structure=="EX") {rho.true=0.2}
                  if(structure=="AR1") {rho.true=0.4}

		for(kk in 1:nsims){
			data=data_generator(n,J,N,Ntest,Kappa,sigma,rho.true,structure)
			YY=data$YY
			Mtest=data$Mtest
			Y=matrix(apply(YY,2,mean),N,1)

	 		#############################################
	 		###		start estimation mhat		###
	 		#############################################
			beta.est=beta%*%Y
			mhat.data=XB1%*%beta.est
			mhat.test=XBtest1%*%beta.est         

			#############################################
			###	start estimation Ghat (jj1,jj2)	###
			#############################################
			Gest=Ghat(YY,mhat.data,Beta,xknots,xtestknots,Ntest)
			phiavgY.data=Gest$phiavgY
			lamdaK.data=Gest$lamdaK
			phiavg.data=Gest$phiavg
                  
			if(select==TRUE){
				Mtest<-mhat.test
				mhat.fix<- mhat.data
              		for(str in 1:length(struct)){
					ST=struct[str]
					rho.est=rho_est(YY,ST) 
					var.zeta=var_zeta(phiavg.data,lamdaK.data,rho.est,J,ST)
					Q.est=generateQ(J,rho.est,ST,phiavg.data,lamdaK.data,var.zeta,alpha)
					Q=Q.est$Q
					#covM=Q.est$covM

					###################################
					# bootstrapping generate new data #
					###################################
	                        YYnew=array(0, c(J,N,n))	
					for(new in 1:bootnum){                    
						#kxinew=generatekxi(covM,lamdaK.data)
						boot.index=sample(1:n,n,replace=T)
						YYnew=YY[,,boot.index]
						#for(i in 1:n){
						#	YYnew[,,i]=matrix(rep(mhat.fix,J),J,dim(mhat.fix)[1],byrow=T)+kxinew[i,,]%*%t(phiavgY.data)
              				#}  
                        		Ynew=matrix(apply(YYnew,2,mean),N,1)
						beta.est=beta%*%Ynew
						mhat=XB1%*%beta.est
						mhat.test=XBtest1%*%beta.est                                
						Gest=Ghat(YYnew,mhat,Beta,xknots,xtestknots,Ntest)
						lamdaK=Gest$lamdaK
						phiavg=Gest$phiavg

						rho.est=rho_est(YYnew,ST)
						var.zeta=var_zeta(phiavg,lamdaK,rho.est,J,ST)
						cov.rate=cov_rate(mhat.test,Mtest,var.zeta,Q,ST,alpha,n)
						conboot[new,,str]=cov.rate
                    		}
                              #cat(colMeans(conboot[ , , str]),"\n","structure=",str,"\n")
					choose_str[str,]=abs(colMeans(conboot[ , , str])-(1-alpha)) 
				} 
                          aa=which(rowMeans(choose_str)==min(rowMeans(choose_str))) 
                        index[kkk,kk]=min(aa)
                         
            
	                 
			}
			if(select==FALSE){
			   for(str in 1:length(struct)){
					ST=struct[str]
					rho.est=rho_est(YY,ST)
                              #print(rho.est) 
					var.zeta=var_zeta(phiavg.data,lamdaK.data,rho.est,J,ST)
					Q.est=generateQ(J,rho.est,ST,phiavg.data,lamdaK.data,var.zeta,alpha)
					Q=Q.est$Q
					covM=Q.est$covM
					covrate[str,kk,]=cov_rate(mhat.test,Mtest,var.zeta,Q,ST,alpha,n)
				   }
			}
              
		}
            cat("true structure", structure, "\n")
           if(select==T){cat("n=", n,"\n")
                  for(str in 1:length(struct)){cat( struct[str], mean(index[kkk,]==str),"\n") } }
          if(select==F){ cat("n=", n,"\n")
                      for(str in 1:length(struct))  cat( struct[str], colMeans(covrate[str,,]),"\n") }

        
	} 

}










