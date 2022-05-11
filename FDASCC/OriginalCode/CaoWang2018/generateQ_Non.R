generateQNon<-function(J,covM,phiavg,lamdaK,var.zeta,alpha,indexka,kappa.select, num){
	
   
	#################
	## generate Q  ##
	#################

	boot=1000
	Xi=rep(0,boot)

	num.select=min(num)
	kappa.s=min(kappa.select)
	
	if (num.select-kappa.s>0) 
{dim1=num.select-kappa.s+1
 ii=1
	Q=matrix(NA,dim1,2)
	for(kap in kappa.s:num.select)	
	{  znorm=matrix(NA,J,kap)
	for(b in 1:boot){     
		zetab=rep(0,Ntest) 
        		for(ka in 1:kap){
        		 znorm[,ka]<- mvrnorm(1,rep(0,J),covM[,,ka])
        		   for(jj1 in 1:J){
          			zetab=zetab+phiavg[,ka,jj1]*znorm[jj1,ka]*lamdaK[ka,jj1]^0.5
        		   }
        		}
		q.res=(zetab/J)/(var.zeta)^(0.5)
		Xi[b]<- max(abs(q.res))
	}
	Q[ii,]=quantile(Xi,c(1-alpha))
	ii=ii+1
	}
	}
	
	if(num.select-kappa.s==0) 
	  { dim1=1
	  Q=matrix(NA,1,2)
	  for(kap in 1:1)	
	  {  znorm=matrix(NA,J,kap)
	  for(b in 1:boot){     
	    zetab=rep(0,Ntest) 
	    for(ka in 1:kap){
	      znorm[,ka]<- mvrnorm(1,rep(0,J),covM[,,ka])
	      for(jj1 in 1:J){
	        zetab=zetab+phiavg[,ka,jj1]*znorm[jj1,ka]*lamdaK[ka,jj1]^0.5
	      }
	    }
	    q.res=(zetab/J)/(var.zeta)^(0.5)
	    Xi[b]<- max(abs(q.res))
	  }
	  Q[1,]=quantile(Xi,c(1-alpha))
	 
	  }
	}
	Qmax=apply(Q, MARGIN = 2, function(x) max(x, na.rm=TRUE))
	list(Q=Qmax)
}

