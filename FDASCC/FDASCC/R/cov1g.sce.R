#' @import gtools
#' @importFrom gtools combinations
#'

cov1g.sce = function(n, N, kappa, Y, muhat, Ghat, GhatY, lamdaK, boot, alpha.grid = c(0.1, 0.05, 0.01)){
  #Counstruct simultaneous confidnece envelope for the covariance funciton in one sample case

  NN = 200  #test points
  nalpha = length(alpha.grid)
  Xi <- rep(0,boot)
  ##eigenfunctions
  phi <- matrix((eigen(Ghat/N, symmetric=T)$vectors[,1:kappa])*sqrt(N),N,kappa)
  phiY <- matrix((eigen(GhatY/NN, symmetric=T)$vectors[,1:kappa])*sqrt(NN),NN,kappa)

  ##get Vhat
  epi <- ((Y-matrix(rep(muhat,n),n,N,byrow=T))%*%sweep(phi,2,(lamdaK)^(-0.5),'*')*N^(-1))^4
  V.hat <- GhatY^2+diag(GhatY)%*%t(diag(GhatY))+(sweep(phiY^2, MARGIN=2, lamdaK, '*'))%*%t(sweep(sweep(phiY^2, MARGIN=2, lamdaK, '*'), MARGIN=2, (colMeans(epi)-3),'*'))

  ### obtain qunatile: resampling ####
  com <- combinations(kappa,2,seq(1:kappa))
  comb <- dim(com)[1]
  Z2 <- matrix(rnorm(kappa*(kappa-1)/2*boot,0,1),boot,kappa*(kappa-1)/2)
  Z1 <- matrix(rnorm(kappa*boot,0,1),boot,kappa)
  for(b in 1:boot){
    zeta <- 0
    for(i in 1:comb){
      zeta <- zeta+((phiY[,(com[i,1])]*lamdaK[(com[i,1])]^0.5)%*%t(phiY[,(com[i,2])]*lamdaK[(com[i,2])]^0.5)+
                      (phiY[,(com[i,2])]*lamdaK[(com[i,2])]^0.5)%*%t(phiY[,(com[i,1])]*lamdaK[(com[i,1])]^0.5))*Z2[b,(com[i,1])]
    }
    for(i in 1:kappa) {
      zeta <- zeta+(phiY[,i]*lamdaK[i]^0.5)%*%t(phiY[,i]*lamdaK[i]^0.5)*((colMeans(epi)[i]-1))^0.5*Z1[b,i]
    }
    Xi[b] <- max(abs(zeta*(V.hat)^(-0.5)))
  }
  Q <- quantile(Xi,1-alpha.grid)#quantiles for three levels

  ###   construct the envelope     ###
  ####################################
  sce.band <- array(NA,dim=c(NN,NN,2,nalpha))
  for(j in 1:nalpha){
    sce.band[,,1,j] <- GhatY-Q[j]*n^(-1/2)*(V.hat)^0.5#lower bound
    sce.band[,,2,j] <- GhatY+Q[j]*n^(-1/2)*(V.hat)^0.5#upper bound
  }
  list(sce = sce.band)
}
