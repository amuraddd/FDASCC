GhatV<-function(YY,Beta,phiavg,Gavg.diag, xtestknots.v,xknots.v,Ntest,lamdaK){
	J=dim(YY)[1]
	N=dim(YY)[2]
	n=dim(YY)[3]
	
	Cij=matrix(Gavg.diag/J,N^2,1,byrow=T)
	workC=Cij[-c((0:(N-1))*(N+1)+1),]
	Beta%*%(workC)->lamda.hat
	xtestknots.v%*%lamda.hat->BB.v
      matrix(BB.v,Ntest,Ntest,byrow=T)->Ghat_testV.diag
	#BBY.v=xknots.v%*%lamda.hat
	#matrix(BBY.v,N,N,byrow=T)->GhatYV.diag

	###############################
	#   fix kappa number ?
	###############################
   
	#phiavgYv=t(phiavgY)%*%(GhatYV.diag/N)/lamdaKY
	phiavgv=t(phiavg)%*%(Ghat_testV.diag/Ntest)/lamdaK
	
	return( phiavgv=phiavgv)
}
    

