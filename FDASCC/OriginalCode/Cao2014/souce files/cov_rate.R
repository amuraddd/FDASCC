cov_rate<- function(mhat.v.test,Mtest.v,V.zeta,Q,structure,alpha,n){
	size=length(alpha)
	con=rep(0,size)

	for(i in 1:size){
		uptest=mhat.v.test+Q[i]*(V.zeta)^0.5/n^0.5
		lowtest=mhat.v.test-Q[i]*(V.zeta)^0.5/n^0.5
   		con[i]=(sum(( Mtest.v>=lowtest)*( Mtest.v<=uptest))==(Ntest))
	}
  
	return(con)
}