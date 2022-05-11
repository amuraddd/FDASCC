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
source("xknots.R")
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
library(refund)
Order1=Order2=3
alpha=c(0.05,0.01)
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
YY1=data.full[c(1:kn),,]
Kappa=kappa=6#2
J<-kn
r1<-r[,c(1:kn)]
n1=dim(YY1)[3]





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
		  phiavgY.data=Gest$phiavgY
		  lamdaK.data=Gest$lamdaK
		  phiavg.data=Gest$phiavg
		  covM=Gest$covM 
		  GMean=Gest$GMean 
		  indexka=Gest$indexka
		  kappa.select=6 #Gest$kappa.select
		  num=Gest$num
		  Q.est=generateQNon(J,covM,phiavg.data,lamdaK.data,GMean,alpha,indexka, kappa.select, num)
		  Q=Q.est$Q
		  
		  mhat.data1+Q[1]*(GMean)^0.5/(n1)^0.5->uptest1
		  mhat.data1-Q[1]*(GMean)^0.5/(n1)^0.5->lowtest1
		  
		  mhat.data1+Q[2]*(GMean)^0.5/(n1)^0.5->uptest2
		  mhat.data1-Q[2]*(GMean)^0.5/(n1)^0.5->lowtest2
		  
		  
			plot(seq(0,(Ntest-1),by=1), (mhat.data1),ylim=c(min(lowtest1-0.05),max(uptest1)+0.05), 
			     col=1,type="l", lty=1, lwd=4,xlab="CCA tract location", ylab=expression(hat(mu)),cex.axis=1.5,cex.lab=1.5, mgp=c(2.5,1,0))
			lines(seq(0,(Ntest-1),by=1), lowtest2, lty=2,col=2,lwd=4)
			lines(seq(0,(Ntest-1),by=1), uptest2, lty=2,col=2,lwd=4)
			
			plot(seq(0,(Ntest-1),by=1), YY1[1,,1],ylim=c(0.25,0.75), 
			     col=1,type="n", lty=1, lwd=1,xlab="CCA tract location", ylab="FA",cex.axis=1.5,cex.lab=1.5, mgp=c(2.5,1,0))
			for (i in 1:kn)
			{ 
			{ for (k in 1:50) lines(seq(0,(Ntest-1),by=1), YY1[i,,k], lty=1,col=1,lwd=1)}}
			
      






