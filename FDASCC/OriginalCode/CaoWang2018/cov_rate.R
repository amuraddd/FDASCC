cov_rate<- function(mhat.test,Mtest,GMean,Q,alpha,n){
	size=length(alpha)
	con=rep(0,size)

	for(i in 1:size){
		uptest=mhat.test+Q[i]*(GMean)^0.5/n^0.5
		lowtest=mhat.test-Q[i]*(GMean)^0.5/n^0.5
   		con[i]=(sum(( Mtest>=lowtest)*( Mtest<=uptest))==(Ntest))
	}
	plot(seq(1/Ntest,1,by=1/Ntest), Mtest, type="l",ylim=c(min(lowtest),max(uptest)))
	lines(seq(1/Ntest,1,by=1/Ntest), lowtest, type="l",col=2)
	lines(seq(1/Ntest,1,by=1/Ntest), uptest, type="l",col=2)
	
	
	return(con)
}