generateQNon2sample<-function(J,covM1,phiavg.data1,lamdaK.data1,GMean1, covM2,phiavg.data2,lamdaK.data2, GMean2,alpha){
	
      kappa=dim(lamdaK.data1)[1]
      znorm1=znorm2=matrix(NA,J,kappa)
	#################
	## generate Q  ##
	#################
      boot=1000
      Xi=rep(0,boot)
      for(b in 1:boot){     
        zetab=rep(0,N) 
      for(ka in 1:kappa){
           znorm1[,ka]<- mvrnorm(1,rep(0,J),covM1[,,ka])
           znorm2[,ka]<- mvrnorm(1,rep(0,J),covM2[,,ka])
        for(jj1 in 1:J){
                   zetab=zetab+n2^0.5*phiavg.data1[,ka,jj1]*znorm1[jj1,ka]*lamdaK.data1[ka,jj1]^0.5+n1^0.5*phiavg.data2[,ka,jj1]*znorm2[jj1,ka]*lamdaK.data2[ka,jj1]^0.5
        }
      }
        q.res=(zetab/J)/(n2*GMean1+n1*GMean2)^(0.5)
        Xi[b]<- max(abs(q.res))
      }
      Q=quantile(Xi,c(1-alpha))
      list(Q=Q)
}