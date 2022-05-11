rho_est<-function(YY,structure,r,mhat){
  J=dim(YY)[1]
  N=dim(YY)[2]
  n=dim(YY)[3]
	A=matrix(0,J,J)
	for(j1 in 1:J){ in1=matrix(rep(r[,j1],N),N,n,byrow=T)
		for(j2 in 1:J){ in2=matrix(rep(r[,j2],N),N,n,byrow=T)
			#A[j1,j2]<- mean(cov(t(YY[j1,,]),t(YY[j2,,])))
			A[j1,j2]<-mean((YY[j1,,]-matrix(rep(mhat,n),N,n,byrow=F)*in1)%*%(t(YY[j2,,])-matrix(rep(mhat,n),n,N,byrow=T)*t(in2))/sum(r[,j1]*r[,j2]))
		}
	}

	if(structure=="AR1"){
		#Xvector= matrix(A[,-J],J*(J-1),1, byrow=T)
		#Yvector= matrix(A[upper.tri(A,diag=T)][,-1],J*(J-1),1, byrow=T)
            B=A[upper.tri(A,diag=T)]
            Xvector= matrix(B[seq(1,J*(J-1)/2)],J*(J-1)/2,1)
		Yvector= matrix(A[upper.tri(A,diag=F)],J*(J-1)/2,1)

		rho.est=(solve( t(Xvector)%*% Xvector)%*% t(Xvector)%*% Yvector)[1,1]
      } 
	if(structure=="EX"){
		rho.est=(sum(A)-sum(diag(A)))/(J*(J-1))/mean(diag(A)) 
	}
	if(structure=="BAND"){ 
		rho.est=mean(diag(A[,-1][-J,]))/mean(diag(A))  
	}
      if(structure=="TOEP"){
		rho.est=rep(0,J-1)
            for (j in 1:(J-2)){
			rho.est[j]=mean(diag(A[-(1:j),][,-((J-(j-1)):J)]))/mean(diag(A))
		}  
		rho.est[J-1]=A[J,1]/mean(diag(A))
	}
	if(structure=="IND"){rho.est=0}
  
	return(rho.est)
}