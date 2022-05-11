Ghat<-function(YY,mhat,Beta,xknots,xtestknots,Ntest){
	J=dim(YY)[1]
	N=dim(YY)[2]
	n=dim(YY)[3]
	Gavg.diag=matrix(0,N,N)
      for(jj1 in 1:J){
		Gavg.diag<- Gavg.diag+(YY[jj1,,]-matrix(rep(mhat,n),N,n,byrow=F))%*%(t(YY[jj1,,])-matrix(rep(mhat,n),n,N,byrow=T))/n
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
	kappa=sum(cumsum(evalues)/sum(evalues)<0.95)+1
	phi1Y<-as.numeric(resY$vectors[,1:kappa])*sqrt(N)
	phi1Y<-matrix(phi1Y,N,kappa)
	phisign=apply(sign(phi1Y)==sign(phi1Y),2,mean)
	phiavgY=t(t(phi1Y)*(2*(phisign>0.5)-1))

      phi1<-as.numeric(res$vectors[,1:kappa])*sqrt(Ntest)
	phi1<-matrix(phi1,Ntest,kappa)
	phisign=apply(sign(phi1Y)==sign(phi1Y),2,mean)
	phiavg=t(t(phi1)*(2*(phisign>0.5)-1))
    
	lamdaK<- evalues[1:kappa]
	lamdaKY<- evaluesY[1:kappa]
     
	list(phiavgY=phiavgY, lamdaKY=lamdaKY,phiavg=phiavg,lamdaK=lamdaK,Gavg.diag=Gavg.diag)
}
    

