library(MASS)
library(splines)

##Order=Order1=Order2=3 
#0: constant 
#1 linear 
#2 quadratic 
#3 cubic
##if using bs, Order must >2

p=4 ## cubic spline 
v=1

set.seed(500)

sigma=0.1
boot=5000


num=c(30,30,30,30, 50,50,50,50,100,100,100,100)
Num=c(50,100,150,200,50,100,150,200,100,200,300,400)



Ns2<-6

result=matrix(0,2,length(num))
alpha1=0.05
alpha2=0.01

Ntest=200
test= seq(1/Ntest,1,1/Ntest)
 
nsims=500 ### no. of replications for coverage frequence

select=rep(0,nsims)
select.knot=rep(0,length(num))

for (jj in 1:length(num))
{
 n<- num[jj]
 N<-Num[jj] 
 
 ns1=seq(floor(n^(1/6)),min(10*floor((N/log(n))^(1/3)),floor(n/4),20),1)
 X=seq(1/N,1,1/N)

#######################################
# True mu and True V
#######################################
 Kappa=8
 phi=phi.v=matrix(0,N,Kappa)
 phi.test=matrix(0,Ntest,Kappa)
  
 for(k in 1:Kappa)
 {
  if(k<5){
  phi[,k]=2^(1/2)*sin(pi*k*X)
  phi.test[,k]=2^(1/2)*sin(pi*k*test)
  }
  if(k>4){
  phi[,k]=2^(1/2)*cos(pi*(k-4)*X)
  phi.test[,k]=2^(1/2)*cos(pi*(k-4)*test)
  }
 }
 for(k in 1:Kappa)
 {
  if(k<5){
  phi.v[,k]=2^(1/2)*pi*k*cos(pi*k*X)
 }
 
 if(k>4){
 phi.v[,k]=-2^(1/2)*pi*k*sin(pi*k*X)
 }
   }
 lamdaK=1/(2^seq(0,(Kappa-1)))
 TrueM=4*X+(2*pi)^(-0.5)/0.1*exp(-(X-0.5)^2/2/0.01)
 TrueM.test=4*test+(2*pi)^(-0.5)/0.1*exp(-(test-0.5)^2/2/0.01)
   
 TrueM.v=4-(X-0.5)/0.01/0.1*(2*pi)^(-0.5)*exp(-(X-0.5)^2/2/0.01)
 TrueM.v.test=4-(test-0.5)/0.01/0.1*(2*pi)^(-0.5)*exp(-(test-0.5)^2/2/0.01)        ##############*************

 #psi=phi%*%diag(sqrt(lamdaK))
 #G=psi%*%t(psi)

 #phi.v=matrix(0,N,Kappa)
 #for(k in 1:Kappa)
 #{
 # phi.v[,k]=2^(1/2)*k*pi*cos(pi*k*X)
 #}
 #psi.v=phi.v%*%diag(sqrt(lamdaK))
 #V=psi.v%*%t(psi.v)
 #diag(V)


###############################
#	Knots for covariance	#
###############################

 knots2=seq(0,1,length=Ns2+2)
 XB1.p2<-bs(X,knots=knots2[-c(1,(Ns2+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))]) 
 XB1.p2.test<-bs(test,knots=knots2[-c(1,(Ns2+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))])
   
 XB.1.p2<-bs(X,knots=knots2[-c(1,(Ns2+2))],degree=p-2,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))])
 XB.1.p2.test<-bs(test,knots=knots2[-c(1,(Ns2+2))],degree=p-2,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))])

 D=diag(-1,Ns2+p)
 D[row(D) == col(D) + 1]<- 1
 D[,-(Ns2+p)]->D2 

 s=1
 t(D2)->A2
 for (i in 1:(p-s-1))
 {
  A2[i,]*(p-s)/(p-(p-s-i+1))->A2[i,]
  A2[(Ns2+p-i-s+1),]*(p-s)/(p-(p-s-i+1))->A2[(Ns2+p-s),]
 }
 t(A2)->D2

 xknots=matrix(0,(N^2),((Ns2+p)*(Ns2+p)))
 xknots.v=matrix(0,(N^2),((Ns2+p)*(Ns2+p)))
 xknots.v.test=matrix(0,(N*Ntest),((Ns2+p)*(Ns2+p)))

 (Ns2+1)*XB.1.p2%*%t(D2)->XB.1.v
 (Ns2+1)*XB.1.p2.test%*%t(D2)->XB.1.v.test

 for(j in 1:N){
  for(k in 1:N){
  	xknots.v[((j-1)*N+k),]<-kronecker(XB1.p2[j,],XB.1.v[k,])
    }
   }

 for(j in 1:N){
  for(k in 1:N){
  	xknots[((j-1)*N+k),]<-kronecker(XB1.p2[k,],XB1.p2[j,])
    }
   }

 solve(t(xknots)%*%(xknots))%*%t(xknots)->Beta

 for(j in 1:N){
  for(k in 1:(Ntest)){
  	xknots.v.test[((j-1)*(Ntest)+k),]<-kronecker(XB1.p2[j,],XB.1.v.test[k,])
    }
   }
 
##############################
#    repetition 500 times    #
##############################

 Ghat=Vhat=matrix(0,N,N)
 C=matrix(0,N,N)
 Cij<-matrix(0,Ntest^2,1,byrow=T)
 #Cij[-c((0:(N-1))*(N+1)+1),]->workCij
 YY=matrix(0,n,N)
 Y=rep(0,N)
 GCV=rep(0,length(ns1))
 con1=con2=rep(0,nsims) 

###########################
#    Start generating Y   #
###########################

 for(kk in 1:nsims)
 {
  kxi=matrix(0,n,Kappa)
  for(j in 1:Kappa)
  {
   kxi[,j]<-rnorm(n,mean=0,sd=sqrt(lamdaK[j]))
  }

  for(j in 1:N)
  { 
  YY[,j]=TrueM[j]
  for(k in 1:Kappa)
  {
   YY[,j]=YY[,j]+kxi[,k]*phi[j,k]
  }
  YY[,j]=YY[,j]+sigma*(rnorm(n))
 }
 matrix(colMeans(YY),N,1)->Y

############################
#   GCV knot selection     #
############################
 
 for(kn in 1:length(ns1))
 {
  Ns1<-ns1[kn]
  knots1=seq(0,1,length=Ns1+2)

############################
#     Spline Basis         #
############################
 
  XB1<-bs(X,knots=knots1[-c(1,(Ns1+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))]) 
  solve(t(XB1)%*%(XB1))%*%t(XB1)->beta  
  XB1%*%beta->Hlamda
  mhat=Hlamda%*%Y
  GCV[kn]=mean((Y-mhat)^2)/(1- sum(diag(Hlamda))/N)^2
 }

 Ns1=ns1[which.min(GCV)]
 Ns1->select[kk]
 #print(Ns1)

############################
#     Use selected knots   #
############################

 knots1=seq(0,1,length=Ns1+2)
 XB1<-bs(X,knots=knots1[-c(1,(Ns1+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))]) 
 XB1.test<-bs(test,knots=knots1[-c(1,(Ns1+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))]) 
 beta=solve(t(XB1)%*%(XB1))%*%t(XB1) 
 bhat=beta%*%Y
 mhat=XB1%*%bhat
 mhat.test=XB1.test%*%bhat

################
 XB.1<-bs(X,knots=knots1[-c(1,(Ns1+2))],degree =p-2,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))])
 XB.1.test<-bs(test,knots=knots1[-c(1,(Ns1+2))],degree=p-2,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))])                     

 D=diag(-1,Ns1+p)
 D[row(D) == col(D) + 1]<- 1
 D[,-(Ns1+p)]->D1 
 t(D1)->A1
  
 s=1
 for (i in 1:(p-s-1))
 {
  A1[i,]*(p-s)/(p-(p-s-i+1))->A1[i,]
  A1[(Ns1+p-i-s+1),]*(p-s)/(p-(p-s-i+1))->A1[(Ns1+p-i-s+1),]
 }
 t(A1)->D1
        
 mhat.v=as.vector((Ns1+1)*XB.1%*%t(D1)%*%bhat)   
 mhat.v.test=as.vector((Ns1+1)*XB.1.test%*%t(D1)%*%bhat)  
    
###############################
#		Get Ghat		#
###############################

 for(j in 1:N)
 {
   for(k in j:N)
   {
     C[j,k]<-C[k,j]<-mean((YY[,j]-mhat[j])*(YY[,k]-mhat[k])) 
   }
 }
 
 Cij<-matrix(C,N^2,1,byrow=T)
 Beta%*%(Cij)->lamda.hat
 xknots.v%*%lamda.hat->BB.v
 xknots.v.test%*%lamda.hat->BB.v.test
  
 matrix(BB.v.test,N,Ntest,byrow=T)->Ghat.v.test
 matrix(BB.v,N,N,byrow=T)->Ghat.v
    
 xknots%*%lamda.hat->BB
 matrix(BB,N,N,byrow=T)->Ghat

#############################
# Get phihat and phi.v.hat  #
#############################

 evalues<-(eigen(Ghat/N,symmetric=T)$value)
 #Kappa=sum(cumsum(evalues)/sum(evalues)<0.95)+1
 Kappa=8
 lamdaK.hat<-evalues[1:Kappa]

 phi.hat<-(eigen(Ghat/N,symmetric=T)$vectors[,1:Kappa])*sqrt(N)
 phi.hat<-matrix(phi.hat,N,Kappa)
 phisign=apply(sign(phi)==sign(phi.hat),2,mean)
 phi.hat=t(t(phi.hat)*(2*(phisign>0.5)-1))
 t(phi.hat)%*%(Ghat.v/N)/lamdaK.hat->phi.v.hat
 t(phi.hat)%*%(Ghat.v.test/N)/lamdaK.hat->phi.v.hat.test

 plot(X,phi.v[,1])
 lines(X,phi.v.hat[1,],type="l",col=3)

############################### 
#		get Vhat		#
###############################
 Vhat=matrix(0,N,N)
 Vhat.test=matrix(0,Ntest,Ntest)
 psi.v.hat=t(phi.v.hat)%*%diag(sqrt(lamdaK.hat))
 psi.v.hat.test=t(phi.v.hat.test)%*%diag(sqrt(lamdaK.hat))
 psi.v=t(phi.v.hat)%*%diag(sqrt(lamdaK.hat))
 Vhat=psi.v.hat%*%t(psi.v.hat)    
 Vhat.test=psi.v.hat.test%*%t(psi.v.hat.test)
 #diag(Vhat.test)

############################### 
#	generate qunatile		# 
###############################

 boot=1000
 zeta=rep(0,boot)
 
 matrix(rnorm(Kappa*boot,0,1),boot,Kappa)->Z1

 Xi=rep(0,boot)
	
 for(b in 1:boot){
  zeta=rep(0,N)
  for(k in 1:Kappa){
   zeta=zeta+phi.v.hat[k,]*lamdaK.hat[k]^0.5*Z1[b,k]
  }
  var.zeta=rep(0,N)
  for(k in 1:Kappa){
   var.zeta=var.zeta+(phi.v.hat[k,]*lamdaK.hat[k]^0.5)^2
  }
  zeta=zeta/(var.zeta)^(0.5)
  Xi[b]<- max(abs(zeta))
 }
 
 quantile(Xi,c(1-alpha1,1-alpha2))->Q
 
 mhat.v.test+Q[1]*n^(-1/2)*(diag(Vhat.test))^0.5->uptest1
 mhat.v.test-Q[1]*n^(-1/2)*(diag(Vhat.test))^0.5->lowtest1

 mhat.v.test+Q[2]*n^(-1/2)*(diag(Vhat.test))^0.5->uptest2
 mhat.v.test-Q[2]*n^(-1/2)*(diag(Vhat.test))^0.5->lowtest2
 
 (sum( ((TrueM.v.test>=lowtest1)*(TrueM.v.test<=uptest1)))== Ntest)->con1[kk]
 (sum( ((TrueM.v.test>=lowtest2)*(TrueM.v.test<=uptest2)))== Ntest)->con2[kk]

 #plot(test,mhat.v.test,type="l", ylim=c(min(lowtest2)-1, max(uptest2)+1), ylab=expression(mu^(1)),xlab="t")
 #lines(test,TrueM.v.test,col=3,lty=4,lwd=2)
 #lines(cbind(test,lowtest1),col=2,lty=2,lwd=2)
 #lines(cbind(test,uptest1),col=2,lty=2,lwd=2)
 #lines(cbind(test,lowtest2),col=4,lty=3,lwd=2)
 #lines(cbind(test,uptest2),col=4,lty=3,lwd=2)

 #print(kk)
 
  }

 mean(con1)->result[1,jj]
 mean(con2)->result[2,jj]
 mean(select)->select.knot[jj]
 
 print(result[,jj])
 }
 #sigma
 sigma
 #num
 num
 #Num
 Num
 #Order1
 result

