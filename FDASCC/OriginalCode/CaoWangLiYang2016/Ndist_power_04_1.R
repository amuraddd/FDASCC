 
rm(list=ls(all=TRUE))
library(MASS)
library(splines)
library(gregmisc)
library(gtools)

n_num <- c(100,200,300,500,800)       # sample size
N_num <- floor(n_num^(0.3)*log(n_num))*4 # points on each curve

A <- c(0.4)
sigma <- 0.01		  # standard deviation of epsilon	
delta <- 0.99             # fraction of selection rate

Order1 <- Order2 <- 1     # spline order for smoothing mean and covariance, linear: 1; cubic: 3
alpha1 <- 0.05
alpha2 <- 0.01

nsims <- 1000                     # no. of replications for coverage frequence
con1 <- con2 <- matrix(0, nsims,length(A))
boot <- 1000                      # resampling number
Xi <- rep(0,boot)                   
NN <- 200                         # test points
test <- seq(1/(NN),1,1/(NN))
Ghat <- V.hat <- matrix(0,NN,NN)
TrueG <- matrix(0,NN,NN)
  
Kappa <- 1000
lamda.true <- rep(0,Kappa)


for (kkk in 1:length(n_num))
{ 
  set.seed(3)
  n <- n_num[kkk]                                 # sample size
  Ns1 <- floor(2*log(n)*n^(1/(2*(Order1+1))))     # knots number (Ns1>>Ns2)
  Ns2 <- floor(4*n^(1/(2*(Order2+1)))*log(log(n)))
  N <- N_num[kkk]                                 # Xij, j <- 1,...,N

  YY <- matrix(0,n,N)
  Y <- rep(0,N)
  X <- seq(1/N,1,1/N)
  C <- matrix(0,N,N)
  Cij <- matrix(0,N^2,1,byrow=T)
  workCij <- Cij[-c((0:(N-1))*(N+1)+1),]

  #################################
  ###      True phi & mu        ###
  #################################
  phi.test <- matrix(0,NN,Kappa)
  phi.true <- matrix(0,N,Kappa)
  for(k in 1: (Kappa/2))
  {
    phi.test[,(2*(k-1)+1)] <- 2^(1/2)*cos(pi*(k)*test)
    phi.test[,(2*k)]       <- 2^(1/2)*sin(pi*(k)*test)
    phi.true[,(2*(k-1)+1)] <- 2^(1/2)*cos(pi*(k)*X)
    phi.true[,(2*k)]       <- 2^(1/2)*sin(pi*(k)*X)
    lamda.true[(2*(k-1)+1)] <- (1/4)^((k-1))
    lamda.true[(2*k)] <- (1/4)^(k)
  }
 
  TrueM <- sin(2*pi*(X-0.5))
  psi <- phi.test%*%diag(sqrt(lamda.true))
  TrueG <- psi%*%t(psi)

  ####################################
  ###     spline basis functions   ###
  ####################################
  XB1 <- bs(X, df=Ns1+Order1+1, degree=Order1, intercept=TRUE)   # B matrix
  XB2 <- bs(X, df=Ns2+Order2+1, degree=Order2, intercept=TRUE)
  XBtest1 <- bs(test, df=Ns1+Order1+1, degree=Order1, intercept=TRUE)
  XBtest2 <- bs(test, df=Ns2+Order2+1, degree=Order2, intercept=TRUE)
  beta <- solve(t(XB1)%*%(XB1))%*%t(XB1)

  xtestknots <- matrix(0,((NN)^2),((Ns2+Order2+1)^2))
  xknots <- matrix(0,(N^2),((Ns2+Order2+1)^2))
  xknots <- kronecker(XB2,XB2)
  xworkknots <- xknots[-c((0:(N-1))*(N+1)+1),]
  Beta <- solve(t(xworkknots)%*%(xworkknots))%*%t(xworkknots)
  xtestknots <- kronecker(XBtest2,XBtest2)
  workknots <- xtestknots[-c((0:(N-1))*(N+1)+1),]
  Beta.test <- solve(t(workknots)%*%(workknots))%*%t(workknots)

  #####################################
  ### 	Find the power 		###
  #####################################
 
  for(j_num in 1:length(A))
  {
   for(kk in 1:nsims)
   {
    kxi <- matrix(0,n,Kappa)
    for(j in 1:Kappa)
    {
      kxi[,j] <- rnorm(n,mean=0,sd=sqrt(lamda.true[j]))
    }
    for(j in 1:N)
    {
      YY[,j]=TrueM[j]
      for(k in 1:Kappa)
      {
        YY[,j] <- YY[,j]+kxi[,k]*phi.true[j,k] 
      }
      YY[,j] <- YY[,j]+sigma*(rnorm(n))+sqrt(A[j_num])*2*cos(3*pi*X[j])*kxi[,5]

    }
    Y <- matrix(colMeans(YY),N,1)
    lamda <- beta%*%Y
    mhat <- XB1%*%lamda
 
    ###################################
    ###          Get Ghat           ###
    ###################################
    C <- tcrossprod(t(YY-matrix(rep(mhat, n), n, N, byrow=T)))/n
    Cij <- matrix(C,N^2,1,byrow=T)
    workCij <- Cij[-c((0:(N-1))*(N+1)+1),]
    lamda.hat <- Beta%*%(workCij)
    BB <- xtestknots%*%lamda.hat
    Ghat <- matrix(BB,NN,NN,byrow=T)
    BBY <- xknots%*%lamda.hat
    GhatY <- matrix(BBY,N,N,byrow=T)

    ####################################
    ### egienvalues & egienfunctions ###
    ####################################
    evalues <- (eigen(GhatY/N,symmetric=T)$value)
    kappa <- sum(cumsum(evalues)/sum(evalues)<delta)+1
    lamdaK <- evalues[1:kappa]
    phiY <- (eigen(GhatY/N, symmetric=T)$vectors[,1:kappa])*sqrt(N)
    phiY <- matrix(phiY,N,kappa)
    phi <- (eigen(Ghat/NN, symmetric=T)$vectors[,1:kappa])*sqrt(NN)
    phi <- matrix(phi,NN,kappa)
    
    ################################
    ##        get Vhat            ##
    ################################
    epi <- matrix(0,n,kappa) #ksihat_ik^4
    for(ii in 1:kappa)
    {
      epi[,ii] <- ((YY-matrix(rep(mhat,n),n,N,byrow=T))%*%phiY[,ii]*(lamdaK[ii])^(-0.5)*N^(-1))^4
    }
    V.hat <- Ghat^2+diag(Ghat)%*%t(diag(Ghat))+(sweep(phi^2, MARGIN=2, lamdaK, '*'))%*%t(sweep(sweep(phi^2, MARGIN=2, lamdaK, '*'), MARGIN=2, (colMeans(epi)-3),'*'))

    ####################################
    ### obtain qunatile: resampling ####
    ####################################
    if(kappa<2)
     {matrix(rnorm(kappa*boot,0,1),boot,kappa)->Z1
      for(b  in 1:boot)
       {
        zeta=(phi[,1]*lamdaK[1]^0.5)%*%t(phi[,1]*lamdaK[1]^0.5)*((colMeans(epi)[1]-1))^0.5*Z1[b,1]
        zeta=zeta*(V.hat)^(-0.5)
        Xi[b]<- max(abs(zeta))
       }
     }
    else{
      com <- combinations(kappa,2,seq(1:kappa))
      comb <- dim(com)[1]
      Z2 <- matrix(rnorm(kappa*(kappa-1)/2*boot,0,1),boot,kappa*(kappa-1)/2)
      Z1 <- matrix(rnorm(kappa*boot,0,1),boot,kappa)
      for(b in 1:boot)
      { 
        zeta <- 0
        for(i in 1:comb)
       {
        zeta <- zeta+((phi[,(com[i,1])]*lamdaK[(com[i,1])]^0.5)%*%t(phi[,(com[i,2])]*lamdaK[(com[i,2])]^0.5)+
                   (phi[,(com[i,2])]*lamdaK[(com[i,2])]^0.5)%*%t(phi[,(com[i,1])]*lamdaK[(com[i,1])]^0.5))*Z2[b,(com[i,1])]
       }
       for(i in 1:kappa)
       {
        zeta <- zeta+(phi[,i]*lamdaK[i]^0.5)%*%t(phi[,i]*lamdaK[i]^0.5)*((colMeans(epi)[i]-1))^0.5*Z1[b,i]
       }
       zeta <- zeta*(V.hat)^(-0.5)
	 Xi[b] <- max(abs(zeta))
    }
    Q <- quantile(Xi,c(1-alpha1,1-alpha2))
   } 
    ####################################
    ###   construct the envelope     ###
    ####################################
    uptest1 <- Ghat+Q[1]*n^(-1/2)*(V.hat)^0.5
    lowtest1 <- Ghat-Q[1]*n^(-1/2)*(V.hat)^0.5
    uptest2 <- Ghat+Q[2]*n^(-1/2)*(V.hat)^0.5
    lowtest2 <- Ghat-Q[2]*n^(-1/2)*(V.hat)^0.5
    con1[kk,j_num] <- (sum((TrueG>=lowtest1)*(TrueG<=uptest1))==(NN)^2)
    con2[kk,j_num] <- (sum((TrueG>=lowtest2)*(TrueG<=uptest2))==(NN)^2)
  }
  cat ("delta=", A[j_num], "\n")
  cat ("n=",n, "\n")
  cat( mean(con1[,j_num]), "&",mean(con2[,j_num]), "&", "\n")
 } 
}


