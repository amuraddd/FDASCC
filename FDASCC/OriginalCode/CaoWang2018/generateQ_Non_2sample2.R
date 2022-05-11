generateQNon2sample<-function(J,covM1,phiavg.data1,lamdaK.data1,GMean1, covM2,phiavg.data2,lamdaK.data2, GMean2,alpha,kappa.select1, num1,kappa.select2, num2){
	
  num.select=min(num1,num2)
  kappa.s=min(kappa.select1,kappa.select2)
  
      znorm1=znorm2=matrix(NA,J,num.select)
	#################
	## generate Q  ##
	#################
      boot=1000
      Xi=rep(0,boot)
  if (num.select-kappa.s>0) 
    {dim1=num.select-kappa.s+1
     ii=1
      Q=matrix(NA,dim1,2)
      for(kap in kappa.s:num.select)	
    {
      for(b in 1:boot){     
        zetab=rep(0,N) 
      for(ka in 1:kap){
           znorm1[,ka]<- mvrnorm(1,rep(0,J),covM1[,,ka])
           znorm2[,ka]<- mvrnorm(1,rep(0,J),covM2[,,ka])
        for(jj1 in 1:J){
                   zetab=zetab+n2^0.5*phiavg.data1[,ka,jj1]*znorm1[jj1,ka]*lamdaK.data1[ka,jj1]^0.5+n1^0.5*phiavg.data2[,ka,jj1]*znorm2[jj1,ka]*lamdaK.data2[ka,jj1]^0.5
        }
      }
        q.res=(zetab/J)/(n2*GMean1+n1*GMean2)^(0.5)
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
    { 
     for(b in 1:boot){     
      zetab=rep(0,N) 
     for(ka in 1:kap){
       znorm1[,ka]<- mvrnorm(1,rep(0,J),covM1[,,ka])
       znorm2[,ka]<- mvrnorm(1,rep(0,J),covM2[,,ka])
    for(jj1 in 1:J){
      zetab=zetab+n2^0.5*phiavg.data1[,ka,jj1]*znorm1[jj1,ka]*lamdaK.data1[ka,jj1]^0.5+n1^0.5*phiavg.data2[,ka,jj1]*znorm2[jj1,ka]*lamdaK.data2[ka,jj1]^0.5
     }
   }
  q.res=(zetab/J)/(n2*GMean1+n1*GMean2)^(0.5)
  Xi[b]<- max(abs(q.res))
  }
  Q[1,]=quantile(Xi,c(1-alpha))
  }
 }

Qmax=apply(Q, MARGIN = 2, function(x) max(x, na.rm=TRUE))
list(Q=Qmax)
}