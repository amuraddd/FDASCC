SCB<-function(YY,alpha1)
 {
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
boot=5000
n=dim(YY)[1]  ## number of sample size
N=dim(YY)[2]  ## number of observation points for each curve
ns1=seq(2,min(floor(n/4),20),1)
X=seq(1/N,1,1/N)
Ghat=Vhat=matrix(0,N,N)
#YY=matrix(0,n,N)     ### YY matrix is the observed data and the format of it is n*N
######################################################################
#    Start generating/reading YY (observed functional data Y(t) )     #
######################################################################
# for(i in 1:n)
#  { YY[i,]=10*sin(X)+(rnorm(N)) }
#######################################################################  
matrix(colMeans(YY),N,1)->Y
GCV=rep(0,length(ns1))
Ns2<- 6     ## number of spline knots used for covariance function estimation
select.knot=rep(0,length(num))

###############################
#	Knots for covariance	#
###############################

 knots2=seq(0,1,length=Ns2+2)
 XB1.p2<-bs(X,knots=knots2[-c(1,(Ns2+2))],degree=p-1,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))])    
 XB.1.p2<-bs(X,knots=knots2[-c(1,(Ns2+2))],degree=p-2,intercept=TRUE,Boundary.knots=knots2[c(1,(Ns2+2))]) 
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
 (Ns2+1)*XB.1.p2%*%t(D2)->XB.1.v

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
 beta=solve(t(XB1)%*%(XB1))%*%t(XB1) 
 bhat=beta%*%Y
 mhat=XB1%*%bhat

################
 XB.1<-bs(X,knots=knots1[-c(1,(Ns1+2))],degree =p-2,intercept=TRUE,Boundary.knots=knots1[c(1,(Ns1+2))])
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
###############################
#		Get Ghat		#
###############################
C=matrix(0,N,N)
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
 matrix(BB.v,N,N,byrow=T)->Ghat.v    
 xknots%*%lamda.hat->BB
 matrix(BB,N,N,byrow=T)->Ghat

#############################
# Get phihat and phi.v.hat  #
#############################

 evalues<-(eigen(Ghat/N,symmetric=T)$value)
 Kappa=sum(cumsum(evalues)/sum(evalues)<0.95)+1
 lamdaK.hat<-evalues[1:Kappa]

 phi.hat<-(eigen(Ghat/N,symmetric=T)$vectors[,1:Kappa])*sqrt(N)
 phi.hat<-t(matrix(phi.hat,N,Kappa))
 #phisign=apply(sign(phi)==sign(phi.hat),2,mean)
 #phi.hat=t(t(phi.hat)*(2*(phisign>0.5)-1))
 phi.v.hat=(phi.hat)%*%(Ghat.v/N)/lamdaK.hat

############################### 
#		get Vhat		#
###############################
 Vhat=matrix(0,N,N)
 Vhat.test=matrix(0,Ntest,Ntest)
 psi.v.hat=t(phi.v.hat)%*%diag(sqrt(lamdaK.hat))
 psi.v=t(phi.v.hat)%*%diag(sqrt(lamdaK.hat))
 Vhat=psi.v.hat%*%t(psi.v.hat)    

###############################################################################
#	generate qunatile for constructing derivative function confidence band	# 
##############################################################################
 boot=1000
 matrix(rnorm(Kappa*boot,0,1),boot,Kappa)->Z1
 Xi=rep(0,boot)
	
 for(b in 1:boot){
  zeta=rep(0,N)
  for(k in 1:Kappa){
   zeta=zeta+phi.v.hat[k,]*lamdaK.hat[k]^0.5*Z1[b,k]
  }
  zeta=zeta/(diag(Vhat))^(0.5)
  Xi[b]<- max(abs(zeta))
 }
  quantile(Xi,c(1-alpha1))->Q1
 #####################################################################
#	generate qunatile for constructing mean function confidence band # 
######################################################################
 matrix(rnorm(Kappa*boot,0,1),boot,Kappa)->Z1
  Xi=rep(0,boot)
  for(b in 1:boot){
  zeta=rep(0,N)
  for(k in 1:Kappa){
   zeta=zeta+phi.hat[k,]*lamdaK.hat[k]^0.5*Z1[b,k]
  }
  zeta=zeta/(diag(Ghat))^(0.5)
  Xi[b]<- max(abs(zeta))
 } 
 quantile(Xi,c(1-alpha1))->Q2
 ########################################################################
 # confidence band for  first order derivative function                 #
 #########################################################################
 mhat.v+Q1*n^(-1/2)*(diag(Vhat))^0.5->uptest1     #### upper bound in 95% confidence band for first order derivative function
 mhat.v-Q1*n^(-1/2)*(diag(Vhat))^0.5->lowtest1    #### lower bound in 95% confidence band for first order derivative function
 ########################################################################
 # confidence band for  mean function                                   #
 #########################################################################
 mhat+Q2*n^(-1/2)*(diag(Ghat))^0.5->uptest2         #### upper bound in 95% confidence band for mean function
 mhat-Q2*n^(-1/2)*(diag(Ghat))^0.5->lowtest2        #### lower bound in 95% confidence band for mean function
 ########################################################################
 # plot of confidence band for  first order derivative function         #
 ######################################################################### 
 plot(X,mhat.v,type="l")
 lines(X,uptest1,type="l",col=3)
 lines(X,lowtest1,type="l",col=3)
  ########################################################################
 # plot of confidence band for  mean function                           #
 #########################################################################
 plot(X,mhat,type="l")
 lines(X,uptest2,type="l",col=2)
 lines(X,lowtest2,type="l",col=2)
 
 list(mhat=mhat,mhat.v=mhat.v,mhatv.up=uptest1,mhat.up=uptest2,mhatv.low=lowtest1,mhat.low=lowtest2)
}