  generateQ<-function(covM1, covM2,var.zeta1,var.zeta2,lamdaK.data1,lamdaK.data2,phiavg.data1,phiavg.data2, n1,n2,alpha,kappa.select1, num1,kappa.select2, num2){
    
    num=min(num1,num2)
    kappa.select=min(kappa.select1,kappa.select2)  
    znorm1=znorm2=matrix(NA,J,num)
    
    boot=1000
    Xi=rep(0,boot)
    ii=1
    
    if(kappa.select<num) {dim1=num-kappa.select+1
    
    Q=matrix(NA,dim1,2)
    for(kap in kappa.select:num)	
    {       
	          for(b in 1:boot){     
	               zetab=rep(0,N) 
	               znorm1<-matrix(mvrnorm(num,rep(0,J),covM1),J, num,byrow=T)
	               znorm2<-matrix(mvrnorm(num,rep(0,J),covM2),J, num,byrow=T)
               for(jj1 in 1:J){
        		            for(ka in 1:kap){
          			       zetab=zetab+n2^0.5*phiavg.data1[,ka]*znorm1[jj1,ka]*lamdaK.data1[ka]^0.5+n1^0.5*phiavg.data2[,ka]*znorm2[jj1,ka]*lamdaK.data2[ka]^0.5
        		                    }
      	                     }
		                   q.res=(zetab/J)/(n2*var.zeta1+n1*var.zeta2)^(0.5)
		                   Xi[b]<- max(abs(q.res))
	                     
    }
    Q[ii,]=quantile(Xi,c(1-alpha))
    ii=ii+1
    }
  }
  
  if(kappa.select==num) {dim1=1
  Q=matrix(NA,1,2)
  for(kap in kappa.select:num)	
  {       
    for(b in 1:boot){     
      zetab=rep(0,N) 

      for(jj1 in 1:J){
        for(ka in 1:kap){
          zetab=zetab+n2^0.5*phiavg.data1[,ka]*znorm1[jj1,ka]*lamdaK.data1[ka]^0.5+n1^0.5*phiavg.data2[,ka]*znorm2[jj1,ka]*lamdaK.data2[ka]^0.5
        }
      }
      q.res=(zetab/J)/(n2*var.zeta1+n1*var.zeta2)^(0.5)
      Xi[b]<- max(abs(q.res))
  }
    Q[1,]=quantile(Xi,c(1-alpha))
  }
  }
  Qmax=apply(Q, MARGIN = 2, function(x) max(x, na.rm=TRUE))
  list(Q=Qmax)
     }