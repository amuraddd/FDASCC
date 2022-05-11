generateQ<-function(J,rho.est,structure,phiavg,lamdaK,var.zeta,alpha, kappa.select, num){
 
     
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
	ii=1
if(kappa.select<num) {dim1=num-kappa.select+1

	   Q=matrix(NA,dim1,2)
    for(kap in kappa.select:num)	
    { 
     	for(b in 1:boot){     
		  zetab=rep(0,Ntest) 
		   znorm<-matrix(mvrnorm(kap,rep(0,J),covM),J,kap,byrow=T)
      	for(jj1 in 1:J){
        		for(ka in 1:kap){
          			zetab=zetab+phiavg[,ka]*znorm[jj1,ka]*lamdaK[ka]^0.5
        		}
      	}
		q.res=(zetab/J)/(var.zeta)^(0.5)
		Xi[b]<- max(abs(q.res))
     }
	  Q[ii,]=quantile(Xi,c(1-alpha))
	  ii=ii+1
    }
 }
	
	if(kappa.select==num) {dim1=1
	  Q=matrix(NA,1,2)
	  for(b in 1:boot){     
	    zetab=rep(0,Ntest) 
	    znorm<-matrix(mvrnorm(kap,rep(0,J),covM),J,kap,byrow=T)
	    for(jj1 in 1:J){
	      for(ka in 1:kap){
	        zetab=zetab+phiavg[,ka]*znorm[jj1,ka]*lamdaK[ka]^0.5
	      }
	    }
	    q.res=(zetab/J)/(var.zeta)^(0.5)
	    Xi[b]<- max(abs(q.res))
	  }
	  Q[1,]=quantile(Xi,c(1-alpha))
	}
	
		Qmax=apply(Q, MARGIN = 2, function(x) max(x, na.rm=TRUE))
	list(Q=Qmax,covM=covM)
}