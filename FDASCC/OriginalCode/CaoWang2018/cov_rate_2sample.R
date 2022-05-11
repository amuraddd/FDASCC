cov_rate2sample<- function(mhat.data1,mhat.data2,var.zeta1,var.zeta2,Q,alpha,n1,n2,N){
	size=length(alpha)
	con=rep(0,size)

	for(i in 1:size){
	  uptest=mhat.data1-mhat.data2+Q[i]*(var.zeta1/n1+var.zeta2/n2)^0.5
	  lowtest=mhat.data1-mhat.data2-Q[i]*(var.zeta1/n1+var.zeta2/n2)^0.5
   		con[i]=(sum(( 0>=lowtest)*(0<=uptest))==N)
	}
#	plot(seq(1/Ntest,1,by=1/Ntest), Mtest, type="l",ylim=c(0,20))
#	lines(seq(1/Ntest,1,by=1/Ntest), lowtest, type="l",col=2)
#	lines(seq(1/Ntest,1,by=1/Ntest), uptest, type="l",col=2)
	
	
	return(con)
}