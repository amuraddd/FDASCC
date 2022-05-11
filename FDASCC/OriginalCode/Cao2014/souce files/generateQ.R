generateQ<-function(J,rho.est,structure,phiavg,lamdaK,var.zeta,alpha){
	
      kappa=length(lamdaK)
      Ntest=dim(phiavg)[1]
	if(structure=="AR1"){
		covM=cor_matrix(J,J,rho.est,"AR1")
      } 
	if(structure=="EX"){
		covM=cor_matrix(J,J,rho.est,"EX")
	}
	if(structure=="BAND"){ 
		covM=cor_matrix(J,J,rho.est,"BAND")
	}
      if(structure=="TOEP"){
		covM=toeplitz(c(1,rho.est))
     }
	if(structure=="IND"){covM=diag(J)}

	#################
	## generate Q  ##
	#################

	boot=1000
	Xi=rep(0,boot)
	for(b in 1:boot){     
		zetab=rep(0,Ntest) 
		znorm<-matrix(mvrnorm(kappa,rep(0,J),covM),J,kappa,byrow=T)
      	for(jj1 in 1:J){
        		for(ka in 1:kappa){
          			zetab=zetab+phiavg[,ka]*znorm[jj1,ka]*lamdaK[ka]^0.5
        		}
      	}
		q.res=(zetab/J)/(var.zeta)^(0.5)
		Xi[b]<- max(abs(q.res))
	}
	Q=quantile(Xi,c(1-alpha))
	list(Q=Q,covM=covM)
}