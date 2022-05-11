Ghat<-function(YY,mhat,Beta,xknots,xtestknots,Ntest,r){
	J=dim(YY)[1]
	N=dim(YY)[2]
	n=dim(YY)[3]
	Gavg.diag=matrix(0,N,N)
      for(jj1 in 1:J){ index=matrix(rep(r[,jj1],N),N,n,byrow=T)
		Gavg.diag<- Gavg.diag+(YY[jj1,,]-matrix(rep(mhat,n),N,n,byrow=F)*index)%*%(t(YY[jj1,,])-matrix(rep(mhat,n),n,N,byrow=T)*t(index))/sum(r[,jj1])
	}

	Cij=matrix(Gavg.diag/J,N^2,1,byrow=T)
	workC=Cij[-c((0:(N-1))*(N+1)+1),]
	Beta%*%(workC)->lamda.hat
	xtestknots%*%lamda.hat->BB
  matrix(BB,Ntest,Ntest,byrow=T)->Ghat_test.diag
	BBY=xknots%*%lamda.hat
	matrix(BBY,N,N,byrow=T)->GhatY.diag

	###############################
	#   fix kappa number ?
	###############################

	resY=eigen(GhatY.diag/N,symmetric=T)
	evaluesY<-(resY$value)

	res=eigen(Ghat_test.diag/Ntest,symmetric=T)
	evalues<-(res$value)
	
	evalues.p=evalues[evalues>0]
	kappa.select=sum(cumsum(evalues.p)/sum(evalues.p)<percent)+1
	num=length(evalues.p[evalues>10^(-5)])
	kappa.num=c(kappa.select,num)
	
	lamdaK<-evalues[1:	num]
	lamdaKY<-evaluesY[1:num]
 
	phi1Y<-(resY$vectors[,1:num])*sqrt(N)
	phi1Y<-matrix(phi1Y,N,	num)
	phisign=apply(sign(phi1Y)==sign(phi1Y),2,mean)
	phiavgY=t(t(phi1Y)*(2*(phisign>0.5)-1))

  phi1<-(res$vectors[,1:	num])*sqrt(Ntest)
	phi1<-matrix(phi1,Ntest,num)
	phisign=apply(sign(phi1Y)==sign(phi1Y),2,mean)
	phiavg=t(t(phi1)*(2*(phisign>0.5)-1))
    
     
	list(phiavgY=phiavgY, lamdaKY=lamdaKY,phiavg=phiavg,lamdaK=lamdaK, kappa.num=kappa.num)
}
    

