##############################
# Using real data "DTI$cca" from "refund" package to construct confidence band. 
# This data set has different repeated visit times (different J). But covM (cov(\xi_{ijk},\xi_{ij'k})) is singular.
# Using structure assumptions
## 2-sample test: gender
## c("TOEP" ,"EX","IND") kn=3, select "TOEP"
## gender =1,2 kn=4, ("EX","IND")  "EX"
##############################

rm(list=ls(all=TRUE))

library(MASS)
library(splines)

### Include source files.
source("corr.R")
source("xknots2.R")
source("Ghat.R")
source("Ghat_Nonpara.R")
source("generateQ_Non_2sample2.R")
source("generateQ.R")
source("cov_rate_2sample.R")
source("cov_rate.R")
source("var_zeta.R")
source("rho_est.R")
########################
## initial setting
########################
library(refund)
Order1=Order2=3
alpha=c(0.05,0.001)
percent=0.95


cca=DTI$cca[complete.cases(DTI$cca), ]
maxid=c(1:length(DTI$ID[complete.cases(DTI$cca)]))[DTI$visit[complete.cases(DTI$cca)]==max(DTI$visit[complete.cases(DTI$cca)])]
nscans.comp=(DTI$Nscans[complete.cases(DTI$cca)])
No8=which(DTI$Nscans[complete.cases(DTI$cca)]==8)  # remove nscan=8, because there is only one nscan=8 sample
#cca=cca[-132,]
nscans.comp[c(No8)]<-length(nscans.comp[c(No8)])
nscans.comp[c(314:319)]<-3
control=which(DTI$case==0)
cca=cca[-control,]

sex=DTI$sex
sex=DTI$sex[complete.cases(DTI$cca) ]
#sex=sex[-132]
sex=sex[-control]

n=length(unique( (DTI$ID[complete.cases(DTI$cca)])[-c(control)]  ))
id=(DTI$ID[complete.cases(DTI$cca)])[-c(control)]  
N=dim(cca)[2]
J=max((DTI$visit[complete.cases(DTI$cca)])[-control] )-1		  ## rep(40,length(n_num))
nscans.comp=nscans.comp[-c(132,control)]
data.full=array(0,c(J,N,n))
r=matrix(0,n,J)
gender=rep(NA,n)
id.track=rep(NA,n)

i=index=1
while (i<=dim(cca)[1]) 
  { Nscan=nscans.comp[i]
    #print(c(i,index))
    id.track[index]=id[i]
       {for(j in 1:Nscan)
           {data.full[j,,index]=cca[i+(j-1),]
            r[index,j]=1}
         gender[index]=sex[i+(Nscan-1)]
        i=i+Nscan
        index=index+1}
  
  }
nm=sum(r)
YY=data.full
Kappa=kappa=2

kn=6
YY1=data.full[c(1:kn),,which(gender==1)] ############# male
Kappa=kappa=2
J<-kn
r1<-r[which(gender==1),c(1:kn)]
n1=dim(YY1)[3]                                ########## n1=66

YY2=data.full[c(1:kn),,which(gender==2)] ############# female
r2<-r[which(gender==2),c(1:kn)]
n2=dim(YY2)[3]                               ############# n2=34



	#############################################
	#####			knots matrix		###
	#############################################
	 
	
	knots=knots2(N,N,order1,order2,n1,J)
	XB1=knots$XB1
	beta=knots$beta
	Beta=knots$Beta
	xknots=knots$xknots
	xtestknots=knots$xtestknots
	XBtest1=knots$XBtest1
	#############################################
	###			generate data		###
	#############################################
             
		  Y1=matrix(apply(YY1,2,mean),N,1)
		  Ntest=N
		  #############################################
		  ###  	start estimation mhat		###
		  #############################################
		  beta.est=beta%*%Y1
		  mhat.data1=XB1%*%beta.est
		  mhat.test1=XBtest1%*%beta.est
		  #############################################
		  ###	start estimation Ghat (jj1,jj2)	###
		  #############################################
		  Gest=GhatNP(YY1,mhat.data1,Beta,xknots,xtestknots,Ntest,r1)
		  phiavgY.data1=Gest$phiavgY
		  lamdaK.data1=Gest$lamdaK
		  phiavg.data1=Gest$phiavg
		  covM1=Gest$covM 
		  GMean1=Gest$GMean 
		  kappa.select1=4 #Gest$kappa.select
		  num1=Gest$num
		  
		  ##########################
		  ###### DATA2#################
		  Y2=matrix(apply(YY2,2,mean),N,1)
		  
		  knots=knots2(N,N,order1,order2,n2,J)
		  XB1=knots$XB1
		  beta=knots$beta
		  Beta=knots$Beta
		  xknots=knots$xknots
		  xtestknots=knots$xtestknots
		  XBtest1=knots$XBtest1
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
		  kappa.select2= 4  #Gest$kappa.select
		  num2=Gest$num
		  
		  Q.est=generateQNon2sample(J,covM1,phiavg.data1,lamdaK.data1,GMean1, covM2,phiavg.data2,lamdaK.data2, GMean2,alpha,kappa.select1, num1,kappa.select2, num2)
		  Q=Q.est$Q
	
			    mhat.data1-mhat.data2+Q[1]*(GMean1/n1+GMean2/n2)^0.5->uptest1
			    mhat.data1-mhat.data2-Q[1]*(GMean1/n1+GMean2/n2)^0.5->lowtest1
			
			    mhat.data1-mhat.data2+Q[2]*(GMean1/n1+GMean2/n2)^0.5->uptest2
			    mhat.data1-mhat.data2-Q[2]*(GMean1/n1+GMean2/n2)^0.5->lowtest2
			
			plot(seq(0,(Ntest-1),by=1), (mhat.data1-mhat.data2),ylim=c(min(lowtest1-0.05),max(uptest1)+0.05), 
			     col=1,type="l", lty=4, lwd=4,xlab="tract location", ylab=expression(hat(mu)[M]-hat(mu)[F]),cex.axis=1.5,cex.lab=1.5,  mgp=c(2.5,1,0))
			lines(seq(0,(Ntest-1),by=1), lowtest1, lty=3, lwd=4,col=2)
			lines(seq(0,(Ntest-1),by=1), uptest1, lty=3, lwd=4,col=2)
			lines(seq(0,(Ntest-1),by=1),rep(0,Ntest),lty=1,lwd=3,col=3)
			


			plot(seq(0,(Ntest-1),by=1), (mhat.data1-mhat.data2),ylim=c(min(lowtest2-0.02),max(uptest2)+0.02), 
			     col=1,type="l", lty=4, lwd=4,xlab="CCA tract location", ylab=expression(hat(mu)[M]-hat(mu)[F]),cex.axis=1.5,cex.lab=1.5,  mgp=c(2.5,1,0))
			lines(seq(0,(Ntest-1),by=1), lowtest2, lty=2, lwd=5,col=2)
			lines(seq(0,(Ntest-1),by=1), uptest2, lty=2, lwd=5,col=2)
			lines(seq(0,(Ntest-1),by=1),rep(0,Ntest),lty=1,lwd=4,col=3)
			
			lines(seq(0,(Ntest-1),by=1), lowtest1, lty=3, lwd=4,col=2)
			lines(seq(0,(Ntest-1),by=1), uptest1, lty=3, lwd=4,col=2)
			lines(seq(0,(Ntest-1),by=1),rep(0,Ntest),lty=1,lwd=3,col=3)

			uptest2




