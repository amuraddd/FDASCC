  generateQ<-function(covM1, covM2,var.zeta1,var.zeta2,lamdaK.data1,lamdaK.data2,phiavg.data1,phiavg.data2, n1,n2,alpha){
	                    
                             boot=1000
	                        Xi=rep(0,boot)
	                        for(b in 1:boot){     
	                  	zetab=rep(0,N) 
		                  znorm1<-matrix(mvrnorm(kappa,rep(0,J),covM1),J,kappa,byrow=T)
                              znorm2<-matrix(mvrnorm(kappa,rep(0,J),covM2),J,kappa,byrow=T)
      	                   for(jj1 in 1:J){
        		            for(ka in 1:kappa){
          			       zetab=zetab+n2^0.5*phiavg.data1[,ka]*znorm1[jj1,ka]*lamdaK.data1[ka]^0.5+n1^0.5*phiavg.data2[,ka]*znorm2[jj1,ka]*lamdaK.data2[ka]^0.5
        		                    }
      	                     }
		                   q.res=(zetab/J)/(n2*var.zeta1+n1*var.zeta2)^(0.5)
		                   Xi[b]<- max(abs(q.res))
	                     }
     	                  Q=quantile(Xi,c(1-alpha))
	          list(Q=Q)
     }