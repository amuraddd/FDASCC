RESCB<-function(YY,alpha)
{ 
rm(list=ls(all=TRUE))
library(MASS)
library(splines)

### Include source files.
source("corr.R")
source("xknots2.R")
source("Ghat.R")
source("GhatV.R")
source("generateQV.R")
source("generateQ.R")
source("new_Q.R")
source("cov_rate.R")
source("V_zeta.R")
source("rho_est.R")
source("Dmatrix1.R")
source("Dmatrix2.R")
source("xknots1.R")

 
  boot=500	### generate supermum
  bootnum=500
  Order1=3
  Order2=3
  struct=c( "IND","EX", "AR1","TOEP")  #c("IND", "EX", "AR1","TOEP")
  conboot=matrix(0,bootnum,length(struct))
  choose_str=matrix(0,length(struct),1)
  index=matrix(0,1,1)
  covrate=array(0,c(length(struct),1,1))
  n=dim(YY)[3]
  N=Ntest=dim(YY)[2]
  J=dim(YY)[1]  ### repeated measurements number
  Y=matrix(apply(YY,2,mean),N,1)
  X=seq(1/N,1,1/N)
	#############################################
	#####	knots matrix and D matrix		###
	#############################################
  D2=Dmatrix2(n,Order2)
	knots=knots2(N,Ntest,order1,order2,n,J,D2)
	Beta=knots$Beta
	xknots=knots$xknots
	xtestknots=knots$xtestknots
  xtestknots.v=knots$xtestknots.v

	#############################################
	###			real data		###
	#############################################
 Ns1= floor(1*log(n)*n^(1/(2*(Order1)))) #GCVknot(Y,ns1,N,Order1) #  
 knots11=knots1(N,Ntest,Ns1,Order1)
 beta=knots11$beta
 XB1=knots11$XB1
 XB.1.test=knots11$XB.1.test
 D1=Dmatrix1(n,Order1,Ns1)
 beta.est=beta%*%Y   
 mhat.data=XB1%*%beta.est
 mhat.v.test=as.vector((Ns1+1)*XB.1.test%*%t(D1)%*%beta.est)
 Mtest.v<-mhat.v.test  
 #############################################
	   ###	start estimation Ghat (jj1,jj2)	###
 #############################################
 Gest=Ghat(YY,mhat.data,Beta,xknots,xtestknots,Ntest)
 phiavg.data=Gest$phiavg
 lamdaK.data=Gest$lamdaK
 Gavg.diag=Gest$Gavg.diag
 GVest=GhatV(YY,Beta,phiavg.data,Gavg.diag, xtestknots.v,xknots.v,Ntest,lamdaK.data)
 phiavgv.data=GVest   
			
 for(kst in 1:length(struct)){
		ST=struct[kst]
		rho.est=rho_est(YY,ST)                          
	###################################
	# bootstrapping generate new data #
	###################################
  YYnew=array(0, c(J,N,n))	
 for(new in 1:bootnum)
	{                    
	boot.index=sample(1:n,n,replace=T)
	YYnew=YY[,,boot.index]
  Ynew=matrix(apply(YYnew,2,mean),N,1)
	beta.est=beta%*%Ynew
	mhat.data=XB1%*%beta.est
	mhat.v.test=as.vector((Ns1+1)*XB.1.test%*%t(D1)%*%beta.est)  
			#############################################
			###	start estimation Ghat (jj1,jj2)	###
			#############################################
	Gest=Ghat(YYnew,mhat.data,Beta,xknots,xtestknots,Ntest)
	phiavg.data=Gest$phiavg
	lamdaK.data=Gest$lamdaK
  Gavg.diag=Gest$Gavg.diag
  GVest=GhatV(YYnew,Beta,phiavg.data,Gavg.diag, xtestknots.v,xknots.v,Ntest,lamdaK.data)
	phiavgv.data=GVest 
	V.zeta=V_zeta(phiavgv.data,lamdaK.data,rho.est,J,ST)
	Q.est=generateQV(J,rho.est,ST,phiavgv.data,lamdaK.data,V.zeta,alpha)
	Q=Q.est$Q
	covM=Q.est$covM
	conboot[new,kst]=cov_rate(mhat.v.test,Mtest.v,V.zeta,Q,ST,alpha,n)
   }					                    
  choose_str[kst,]=abs(mean(conboot[, kst])-(1-alpha)) 
	} 
  aa=which(rowMeans(choose_str)==min(rowMeans(choose_str))) 
  index=min(aa)
	ST=struct[index]
	rho.est=rho_est(YY,ST)                          
	V.zeta=V_zeta(phiavgv.data,lamdaK.data,rho.est,J,ST)
	Q.est=generateQV(J,rho.est,ST,phiavgv.data,lamdaK.data,V.zeta,alpha)
	Q2=Q.est$Q
	covM=Q.est$covM

	var.zeta=V_zeta(t(phiavg.data),lamdaK.data,rho.est,J,ST)
	Q.est1=generateQ(J,rho.est,ST,t(phiavg.data),lamdaK.data,var.zeta,alpha)
	Q1=Q.est1$Q
	
				 mhat.data+Q1*(var.zeta)^0.5/(n)^0.5->uptest1
         mhat.data-Q1*(var.zeta)^0.5/(n)^0.5->lowtest1
         mhat.v.test+Q2*(V.zeta)^0.5/(n)^0.5->uptest2
         mhat.v.test-Q2*(V.zeta)^0.5/(n)^0.5->lowtest2
         
         plot(mhat.data,type="n",xlim=c(0,1), xlab="Age", ylab=expression(paste(hat(mu ))),cex.axis=1,cex.lab=1.1,cex.main=1.2)
         lines(cbind(X,lowtest1),col=3,lty=2,lwd=3)
         lines(cbind(X,uptest1),col=3,lty=2,lwd=3)
         lines(cbind(X,mhat.data),col=1,lty=1,lwd=1)
         
         plot(mhat.v.test,type="n",xlim=c(0,1), xlab="Age", ylab=expression(paste(hat(mu )^(1))),cex.axis=1,cex.lab=1.1,cex.main=1.2)
         lines(cbind(X,lowtest2),col=2,lty=2,lwd=3)
         lines(cbind(X,uptest2),col=2,lty=2,lwd=3)
         lines(cbind(X,mhat.v.test),col=1,lty=1,lwd=1)
  
  list(mhat=mhat.data,mhat.v=mhat.v.test,mhatv.up=uptest2,mhat.up=uptest1,mhatv.low=lowtest2,mhat.low=lowtest1)
}



        