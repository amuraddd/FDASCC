GhatNP<-function(YY,mhat.data,Beta,xknots,xtestknots,Ntest,r){
	J=dim(YY)[1]
	N=dim(YY)[2]
	n=dim(YY)[3]
	Gj=array(NA,c(N,N,J))
	Ghat_test=array(NA,c(Ntest,Ntest,J))
	phiavgY=array(NA,c(N,kappa,J))
	phiavg=array(NA,c(Ntest,kappa,J))
	lamdaK=lamdaKY=matrix(NA,kappa,J)
	kappa.select=num=matrix(NA,1,J)
	epi <- array(NA,c(n,kappa,J)) 
	evalues=evaluesY=matrix(NA,kappa,J)
	covM=array(NA,c(J,J,kappa))
	indexka=matrix(0,kappa,2)
	
  for(jj1 in 1:J)
  { index=matrix(rep(r[,jj1],N),N,n,byrow=T)
		Gj[,,jj1]<-((YY[jj1,,]-matrix(rep(mhat.data,n),N,n,byrow=F) )*index )%*%(t(YY[jj1,,])-matrix(rep(mhat.data,n),n,N,byrow=T)*t(index))/sum(r[,jj1])


		Cij=matrix(Gj[,,jj1],N^2,1,byrow=T)
		workC=Cij[-c((0:(N-1))*(N+1)+1),]
		Beta%*%(workC)->lamda.hat
		xtestknots%*%lamda.hat->BB
		matrix(BB,Ntest,Ntest,byrow=T)->Ghat_test[,,jj1]
		BBY=xknots%*%lamda.hat
		matrix(BBY,N,N,byrow=T)->GhatY
		
		###############################
		#   fix kappa number ?
		###############################
		
		resY=eigen(GhatY/N,symmetric=T)
		evaluesY[,jj1]<-(resY$value)[1:kappa]
		
		res=eigen(Ghat_test[,,jj1]/Ntest,symmetric=T)
		evalues[,jj1]<-(res$value)[1:kappa]
		
		evalues.p=evalues[evalues[,jj1]>0]
		kappa.select[,jj1]=sum(cumsum(evalues.p)/sum(evalues.p)<percent)+1
		num[,jj1]=length(evalues.p[evalues>10^(-5)])
		
		phi1Y<-(resY$vectors[,1:kappa])*sqrt(N)
		phi1Y<-matrix(phi1Y,N,kappa)
		phisign=apply(sign(phi1Y)==sign(phi1Y),2,mean)
		phiavgY[,,jj1]=t(t(phi1Y)*(2*(phisign>0.5)-1))
		
		phi1<-(res$vectors[,1:kappa])*sqrt(Ntest)
		phi1<-matrix(phi1,Ntest,kappa)
		phisign=apply(sign(phi1Y)==sign(phi1Y),2,mean)
		phiavg[,,jj1]=t(t(phi1)*(2*(phisign>0.5)-1))
		
		lamdaK[,jj1]<- evalues[1:kappa,jj1]
		lamdaKY[,jj1]<- evaluesY[1:kappa,jj1]
		
	for(ii in 1:kappa)
	{
	  epi[,ii,jj1] <- ((t(YY[jj1,,])-matrix(rep(mhat.data,n),n,N,byrow=T))*t(index))%*%phiavgY[,ii,jj1]*(lamdaKY[ii,jj1])^(-0.5)*N^(-1)
	}
	
  
  }
# 


 
 Gjj=matrix(0,Ntest,1)

 for (jj1 in 1:J)
 { 
   for (jj2 in 1:J)
   { rnm=sum(r[,jj1]*r[,jj2])
     for (k1 in 1:kappa)
     {
       for(k2 in 1:kappa) 
       {Gjj=Gjj+1/rnm*sum(epi[,k1,jj1]*epi[,k2,jj2])*phiavg[,k1,jj1]*(lamdaK[k1,jj1])^(0.5)*(phiavg[,k2,jj2]*(lamdaK[k2,jj2])^(0.5))/(J^2)}
     }
   }
 }
 
 
 
GMean<-Gjj
 
 
 for (jj1 in 1:J)
 {
   for (jj2 in (jj1):J)
   {rnm=sum(r[,jj1]*r[,jj2])
    {for (k1 in 1:kappa)
      covM[jj1,jj2,k1]=covM[jj2,jj1,k1]= sum((epi[,k1,jj1]-sum(epi[,k1,jj1])/sum(r[,jj1]))*r[,jj1]* (epi[,k1,jj2]-sum(epi[,k1,jj2])/sum(r[,jj2]))*r[,jj2])/(rnm-1)   
      #covM[jj1,jj2,k1]=covM[jj2,jj1,k1]=sum(epi[,k1,jj1]*epi[,k1,jj2])/rnm -sum(epi[,k1,jj1])/sum(r[,jj1])*sum(epi[,k1,jj2])/sum(r[,jj2])
      }
   }
 }

for (k1 in 1:kappa)
 { #M<-covM[,,k1]
   eig=eigen(covM[,,k1],T)
   val=eig$values
   vec=eig$vectors
   num=sum(val>0)
  if(num<J) 
    { valpos=val[1:num]
    covM[,,k1]<- vec[,1:num]%*% diag(valpos)%*%t(vec[,1:num])
    indexka[k1,]=c(k1,num)
    }
 }

list(phiavgY=phiavgY, lamdaKY=lamdaKY,phiavg=phiavg,lamdaK=lamdaK,epi=epi,GMean=GMean,covM=covM,indexka=indexka, 	kappa.select=	kappa.select, num=num)
}


