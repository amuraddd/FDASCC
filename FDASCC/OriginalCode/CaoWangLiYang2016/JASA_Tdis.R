#############################################################
###		            Confidence Envelopes		                ###
### 		kxi_k ~ rescaled t(df=10), variance lamda_k  		  ###
### 		              lamda_k=k^(-5)   	                  ###
#############################################################

rm(list=ls(all=TRUE))

library(MASS)
library(splines)
library(gregmisc)
 
n_num <- c(50,100,200,300,500,800)        # sample size
N_num <- floor(n_num^(0.3)*log(n_num))*4  # points on each curve

set.seed(3)
sigma <- 0.01          # standard deviation of epsilon
delta <- 0.99          # fraction of selection rate

Order1 <- Order2 <- 1  # spline order used for smoothing mean and covariance functions, linear spline: 1; cubic spline: 3
alpha1 <- 0.05
alpha2 <- 0.01
degr <- 10             # t distribution degree of freedom
  
nsims <- 1000                       # no. of replications for coverage frequence
con1 <- con2 <- matrix(0, nsims,length(n_num))
boot <- 1000                        ## resampling number
Xi <- rep(0,boot)                   
NN <- 200                           ## test points
test <- seq(1/(NN),1,1/(NN))
Ghat <- V.hat <- matrix(0,NN,NN)
TrueG <- matrix(0,NN,NN)
  
Kappa <- 1000
lamda.true <- rep(0,Kappa)

for(kkk in 1:length(n_num))
{
  n <- n_num[kkk]                                 ### sample size
  Ns1 <- floor(2*log(n)*n^(1/(2*(Order1+1))))     ## knots number (Ns1>>Ns2)
  Ns2 <- floor(4*n^(1/(2*(Order2+1)))*log(log(n)))
  N <- N_num[kkk]                                 #### Xij, j=1,...,N

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
    lamda.true[(2*(k-1)+1)] <- k^(-5)
    lamda.true[(2*k)] <- k^(-5)
  }
 
  TrueM <- sin(2*pi*(X-0.5))
  psi <- phi.test%*%diag(sqrt(lamda.true))
  TrueG <- psi%*%t(psi)

  ##################################################################################
  ###  spline basis functions
  ##################################################################################
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

  #########################################
  ###      starting the simulation      ###
  #########################################
  for(kk in 1:nsims)
  {
    kxi <- matrix(0,n,Kappa)
    kxi <- matrix(rt(n*Kappa, degr)*sqrt((degr-2)/degr),n,Kappa,byrow=F)* matrix(rep(sqrt(lamda.true),n),n,Kappa,byrow=T)  
    YY <- matrix(rep(TrueM,n),n,N,byrow=T)+kxi %*% t(phi.true)+sigma* rt(n, degr)*sqrt((degr-2)/degr)
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
    V.hat <- Ghat^2+diag(Ghat)%*%t(diag(Ghat))+(sweep(phi^2, MARGIN=2, lamdaK, '*'))%*% t(sweep(sweep(phi^2, MARGIN=2, lamdaK, '*'), MARGIN=2, (colMeans(epi)-3),'*'))
    
    ####################################
    ### obtain qunatile: resampling ####
    ####################################
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
    
   ####################################
   ###    construct  envelopes      ###
   ####################################
   uptest1 <- Ghat+Q[1]*n^(-1/2)*(V.hat)^0.5
   lowtest1 <- Ghat-Q[1]*n^(-1/2)*(V.hat)^0.5
   uptest2 <- Ghat+Q[2]*n^(-1/2)*(V.hat)^0.5
   lowtest2 <- Ghat-Q[2]*n^(-1/2)*(V.hat)^0.5
   con1[kk,kkk] <- (sum((TrueG>=lowtest1)*(TrueG<=uptest1))==(NN)^2)
   con2[kk,kkk] <- (sum((TrueG>=lowtest2)*(TrueG<=uptest2))==(NN)^2)
  }
  cat( mean(con1[,kkk]), "&", mean(con2[,kkk]), "&", "\n")
}

res <- cbind(colMeans(con1), colMeans(con2))
write.table(res, paste("Tdis_", as.character(sigma),"_",as.character(Order1),".txt", sep=""))




