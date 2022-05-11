knots2=function(N,Ntest,order1,order2,n,J){
	Ns1=floor(1*log(n)*n^(1/(2*(Order1+1))))	## knots number (Ns1>>Ns2)
	Ns2=floor(2*n^(1/(2*(Order2+1)))*(log(n))^(0.5))

#################################################################################
#######  knots matrix														#####
#################################################################################
	knots1=seq(0,1,length=Ns1+2)
	knots2=seq(0,1,length=Ns2+2)

	X=seq(1/N,1,1/N)	
	XB1<-bs(X,knots=knots1[-c(1,(Ns1+2))],degree=Order1,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))])
	XB2<-bs(X,knots=knots2[-c(1,(Ns2+2))],degree=Order1,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))]) 

	test=seq(1/(Ntest),1,1/(Ntest))
	XBtest1<-bs(test,knots=knots1[-c(1,(Ns1+2))],degree=Order1,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))])
	XBtest2<-bs(test, knots=knots2[-c(1,(Ns2+2))],degree=Order2,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))]) 

	beta<-solve(t(XB1)%*%(XB1))%*%t(XB1)
 
	xtestknots=matrix(0,((Ntest)^2),((Ns2+Order2+1)^2))
	xknots=matrix(0,(N^2),((Ns2+Order2+1)^2))

	xknots=kronecker(XB2,XB2)
	xworkknots<-xknots[-c((0:(N-1))*(N+1)+1),]

	Beta<-solve(t(xworkknots)%*%(xworkknots))%*%t(xworkknots)

	xtestknots<-kronecker(XBtest2,XBtest2)

	list(beta=beta,Beta=Beta,XB1=XB1,xknots=xknots,xtestknots=xtestknots,XBtest1=XBtest1,XBtest2=XBtest2)
}
