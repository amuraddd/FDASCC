###  Data: YY J*N*n   J: multi level.  N: observation number for each sample. n: sample size
###		  Y  1*N	 average YY
###		  X  1:N	 obervation points	
###		  phi.true   N*J:Kappa	(N*J), Kappa
###		  kxi		 n*J*Kappa	(n*J), Kappa
### 
###		  lamdaK=array(0,c(20,J,J))	## in case kappa is large. set 20 as the largest kappa
###		  phi=array(0,c(Ntest,20,J))
###  Generate data with different J

data_generator=function(n,J,N,Ntest,Kappa,sigma,rho.true,structure,mpr){  #mpr: missing probability
  YY=array(0,c(J,N,n))	### n matrices, each matrix: J*N 
  Y=rep(0,N)
  X=seq(1/N,1,1/N)
  test=seq(1/(Ntest),1,1/(Ntest))
  kxi=array(0,c(n,J,Kappa))	### n matrices, each matrix: J*Kappa
  phi.true=matrix(0,N,Kappa)
  select.matrix=matrix(NA,n,J)
  #######################################
  # True mu
  #######################################
  TrueM=10+sin(2*pi*(X-0.5))
  Mtest=10+sin(2*pi*(test-0.5))
  
  #######################################
  # True phi and True lambda
  #######################################
  lamda.true=rep(0,Kappa)
  for(k in 1: (Kappa/2)){
    phi.true[,(2*(k-1)+1)]=2^(1/2)*cos(pi*k*X)
    phi.true[,(2*k)]		=2^(1/2)*sin(pi*k*X)
    #phi.test[,(2*(k-1)+1)]=2^(1/2)*cos(pi*k*test)
    #phi.test[,(2*k)]		=2^(1/2)*sin(pi*k*test)
    lamda.true[(2*(k-1)+1)]= lamda.true[(2*k)]=(1/2)^(2*k)
  }
  
  ###########################################################
  # True epi, kxi														###
  ###########################################################
  for(ka in 1:Kappa){
    Sigma.kxi=cor_matrix(J,J,rho.true,structure)*lamda.true[ka]
    kxi[,,ka]<-mvrnorm(n,rep(0,J),Sigma.kxi)
  }
  
  Sigma.epi= cor_matrix(J,J,rho.true,structure)*sigma^2
  
  ###########################################################
  # True Y data															###
  ############################################################
  for(i in 1:n){
    epi<-mvrnorm(N,rep(0,J),Sigma.epi)
    YY[,,i]=matrix(rep(TrueM,J),J,N,byrow=T)+kxi[i,,]%*%t(phi.true)+ t(epi) 
    #index=floor(runif(1,floor(J/2),J))
    index=rbinom(J,1,mpr) #rbinom(J,1,mpr)
    select=matrix(rep(index, N), byrow=F, J,N)
    YY[,,i]=YY[,,i]*select
    select.matrix[i,]=index
  }
  
  list(YY=YY,Mtest=Mtest,select=select.matrix)
}