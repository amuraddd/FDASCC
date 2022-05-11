knots2=function(N,Ntest,order1,order2,n,J,D2){
	#Ns1=floor(1*log(n)*n^(1/(2*(Order1+1))))	## knots number (Ns1>>Ns2)
	Ns2=floor(2*n^(1/(2*(Order2+1)))*log(log(n)))

#################################################################################
#######  knots matrix														#####
#################################################################################
	#knots1=seq(0,1,length=Ns1+2)
	knots2=seq(0,1,length=Ns2+2)

	X=seq(1/N,1,1/N)	
	#XB1<-bs(X,knots=knots1[-c(1,(Ns1+2))],degree=Order1,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))])
	XB2<-bs(X,knots=knots2[-c(1,(Ns2+2))],degree=Order1,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))]) 

	test=seq(1/(Ntest),1,1/(Ntest))
	#XBtest1<-bs(test,knots=knots1[-c(1,(Ns1+2))],degree=Order1,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))])
	XBtest2<-bs(test, knots=knots2[-c(1,(Ns2+2))],degree=Order2,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))]) 
 
      #XB.1.test<-bs(test,knots=knots1[-c(1,(Ns1+2))],degree=Order1-1,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))])                     
      XB.2.test<-bs(test,knots=knots2[-c(1,(Ns2+2))],degree=Order2-1,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))])                     
      XB.2<-bs(X,knots=knots2[-c(1,(Ns2+2))],degree=Order2-1,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))])                     
  

	#beta<-solve(t(XB1)%*%(XB1))%*%t(XB1)
 
	xtestknots=matrix(0,((Ntest)^2),((Ns2+Order2+1)^2))
	xknots=matrix(0,(N^2),((Ns2+Order2+1)^2))

	xknots=kronecker(XB2,XB2)
	xworkknots<-xknots[-c((0:(N-1))*(N+1)+1),]

	Beta<-solve(t(xworkknots)%*%(xworkknots))%*%t(xworkknots)

	xtestknots<-kronecker(XBtest2,XBtest2)

      XB.1.v.test<-(Ns2+1)*XB.2.test%*%t(D2)
      XB.1.v<-(Ns2+1)*XB.2%*%t(D2)
      xtestknots.v<-kronecker(XBtest2, XB.1.v.test)
      xknots.v<-kronecker(XB2, XB.1.v)

	list(Beta=Beta,xknots=xknots,xtestknots=xtestknots,XBtest2=XBtest2,  xtestknots.v= xtestknots.v, Ns2=Ns2)
}
