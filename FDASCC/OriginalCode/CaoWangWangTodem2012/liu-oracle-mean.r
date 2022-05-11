library(MASS)
library(splines)

 ##Order=Order1=Order2=3 
 #0: constant 
 #1 linear 
 #2 quadratic 
 #3 cubic
 ##if using bs, Order must >2

 p=4 ## cubic spline 
 v=0

 set.seed(500)

 sigma=0.1
 boot=5000

 Kappa=8
 lamdaK=1/2^seq(0,(Kappa-1))*2
 #n=num=c(100)
 #N=200

 #num=c(30,30,50,50,100,100,200,200,300,300)
 #Num=c(30,60,50,100,100,200,200,400,300,500)

 num=c(30,30,30,30,50,50,50,50,100,100,100,100)
 Num=c(50,100,150,200,50,100,150,200,100,200,300,400)
 
 
 result=matrix(0,2,length(num))
 alpha1=0.05
 alpha2=0.01

 Ntest=200
 test= seq(1/Ntest,1,1/Ntest)
 
 nsims=1000 ### no. of replications for coverage frequence

 select=rep(0,nsims)

 select.knot=rep(0,length(num))

 for (jj in 1:length(num))
 {
 n<- num[jj]
 N<-Num[jj] 
 
 ns1=seq(floor(n^(1/6)),min(10*floor((N/log(n))^(1/3)),floor(n/4),20),1)

 X=seq(1/N,1,1/N)


#######################################
# True mu 
#######################################
 
 TrueM=4*X+(2*pi)^(-0.5)/0.1*exp(-(X-0.5)^2/2/0.01)
 TrueM.test=4*test+(2*pi)^(-0.5)/0.1*exp(-(test-0.5)^2/2/0.01)
   
#######################################
# True G 
#######################################


 phi=matrix(0,N,Kappa)
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

 psi.test=phi.test%*%diag(sqrt(lamdaK))
 G.test=psi.test%*%t(psi.test)
 #diag(G.test)
 

##############################
#    repetition 500 times    #
##############################

  C=matrix(0,N,N)
  Cij<-matrix(0,Ntest^2,1,byrow=T)
  #Cij[-c((0:(N-1))*(N+1)+1),]->workCij
  YY=matrix(0,n,N)
  
  Y=rep(0,N) 

  GCV=rep(0, length(ns1))
  con1=con2=rep(0, nsims)

###########################
#    Start generating Y
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
 
for (kn in 1:length(ns1))
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

############################
#     Use selected knots   #
############################

  knots1=seq(0,1,length=Ns1+2)

  XB1<-bs(X,knots=knots1[-c(1,(Ns1+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))]) 
  XB1.test<-bs(test,knots=knots1[-c(1,(Ns1+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))]) 

  solve(t(XB1)%*%(XB1))%*%t(XB1)->beta  
  beta%*%Y->bhat
  mhat=XB1%*%bhat
  mhat.test=XB1.test%*%bhat
      
############################# 
# generate qunatile 
#############################


 Z1=matrix(rnorm(Kappa*boot,0,1),boot,Kappa)
 Xi=rep(0,boot)
	
 for(b in 1:boot){
  zeta=rep(0,N)
  for(k in 1:Kappa){
   zeta=zeta+phi[,k]*lamdaK[k]^0.5*Z1[b,k]
  }
  var.zeta=rep(0,N)
  for(k in 1:Kappa){
   var.zeta=var.zeta+(phi[,k]*lamdaK[k]^0.5)^2
  }
  zeta=zeta/(var.zeta)^(0.5)
  Xi[b]<- max(abs(zeta))
 }
 
 quantile(Xi,c(1-alpha1,1-alpha2))->Q
 
 mhat.test+Q[1]*n^(-1/2)*(diag(G.test))^0.5->uptest1
 mhat.test-Q[1]*n^(-1/2)*(diag(G.test))^0.5->lowtest1

 mhat.test+Q[2]*n^(-1/2)*(diag(G.test))^0.5->uptest2
 mhat.test-Q[2]*n^(-1/2)*(diag(G.test))^0.5->lowtest2

 
 (sum( ((TrueM.test>=lowtest1)*(TrueM.test<=uptest1)))== Ntest)->con1[kk]
 (sum( ((TrueM.test>=lowtest2)*(TrueM.test<=uptest2)))== Ntest)->con2[kk]

 plot(test,mhat.test,type="l", ylim=c(min(lowtest2)-1, max(uptest2)+1), ylab=expression(mu),xlab="t")
 lines(test,TrueM.test,col=3,lty=4,lwd=2)
 lines(cbind(test,lowtest1),col=2,lty=2,lwd=2)
 lines(cbind(test,uptest1),col=2,lty=2,lwd=2)
 lines(cbind(test,lowtest2),col=4,lty=3,lwd=2)
 lines(cbind(test,uptest2),col=4,lty=3,lwd=2)

 # print(kk)
 
  }

 mean(con1)->result[1,jj]
 mean(con2)->result[2,jj]
 mean(select)->select.knot[jj]
 
 print(result[,jj])
 }

 sigma

 num
 Num

select.knot
 
 #Order1
 result

