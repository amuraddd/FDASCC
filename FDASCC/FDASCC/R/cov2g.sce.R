#' @import gtools
#' @importFrom gtools combinations
#'

cov2g.sce = function(na, nb, N, kappa.a, kappa.b, Ya, Yb, muhat.a, muhat.b, Ghat.a, Ghat.b,  GhatY.a, GhatY.b, lamdaK.a, lamdaK.b, rab, boot, alpha.grid = c(0.1, 0.05, 0.01)){
  #Counstruct simultaneous confidnece envelope for the difference of two covariance funcitons


  NN = 200  #test points
  nalpha = length(alpha.grid)
  Xi <- rep(0,boot)
  ##eigenfunctions

  ##group1
  phi.a <- matrix((eigen(Ghat.a/N, symmetric=T)$vectors[,1:kappa.a])*sqrt(N),N,kappa.a)
  phiY.a <- matrix((eigen(GhatY.a/NN, symmetric=T)$vectors[,1:kappa.a])*sqrt(NN),NN,kappa.a)

  ##group2
  phi.b <- matrix((eigen(Ghat.b/N, symmetric=T)$vectors[,1:kappa.b])*sqrt(N),N,kappa.b)
  phiY.b <- matrix((eigen(GhatY.b/NN, symmetric=T)$vectors[,1:kappa.b])*sqrt(NN),NN,kappa.b)

  ##get Vhat

  ##group1

  epi.a <- ((Ya-matrix(rep(muhat.a,na),na,N,byrow=T))%*%sweep(phi.a,2,(lamdaK.a)^(-0.5),'*')*N^(-1))^4
  V.hat.a <- GhatY.a^2+diag(GhatY.a)%*%t(diag(GhatY.a))+(sweep(phiY.a^2, MARGIN=2, lamdaK.a, '*'))%*%t(sweep(sweep(phiY.a^2, MARGIN=2, lamdaK.a, '*'), MARGIN=2, (colMeans(epi.a)-3),'*'))

  ##group2
  epi.b <- ((Yb-matrix(rep(muhat.b,nb),nb,N,byrow=T))%*%sweep(phi.b,2,(lamdaK.b)^(-0.5),'*')*N^(-1))^4
  V.hat.b <- GhatY.b^2+diag(GhatY.b)%*%t(diag(GhatY.b))+(sweep(phiY.b^2, MARGIN=2, lamdaK.b, '*'))%*%t(sweep(sweep(phiY.b^2, MARGIN=2, lamdaK.b, '*'), MARGIN=2, (colMeans(epi.b)-3),'*'))



  ### obtain qunatile: resampling ####
  com.a <- combinations(kappa.a,2,seq(1:kappa.a))
  com.b <- combinations(kappa.b,2,seq(1:kappa.b))
  comb.a <- dim(com.a)[1]
  comb.b <- dim(com.b)[1]
  Z2.a <- matrix(rnorm(kappa.a*(kappa.a-1)/2*boot,0,1),boot,kappa.a*(kappa.a-1)/2)
  Z2.b <- matrix(rnorm(kappa.b*(kappa.b-1)/2*boot,0,1),boot,kappa.b*(kappa.b-1)/2)
  Z1.a <- matrix(rnorm(kappa.a*boot,0,1),boot,kappa.a)
  Z1.b <- matrix(rnorm(kappa.b*boot,0,1),boot,kappa.b)
  for(b in 1:boot){
    zeta.a <- 0
    zeta.b <- 0
    for(i in 1:comb.a){
      zeta.a <- zeta.a+((phiY.a[,(com.a[i,1])]*lamdaK.a[(com.a[i,1])]^0.5)%*%t(phiY.a[,(com.a[i,2])]*lamdaK.a[(com.a[i,2])]^0.5)+
                      (phiY.a[,(com.a[i,2])]*lamdaK.a[(com.a[i,2])]^0.5)%*%t(phiY.a[,(com.a[i,1])]*lamdaK.a[(com.a[i,1])]^0.5))*Z2.a[b,(com.a[i,1])]
    }
    for(i in 1:comb.b){
      zeta.b <- zeta.b+((phiY.b[,(com.b[i,1])]*lamdaK.b[(com.b[i,1])]^0.5)%*%t(phiY.b[,(com.b[i,2])]*lamdaK.b[(com.b[i,2])]^0.5)+
                          (phiY.b[,(com.b[i,2])]*lamdaK.b[(com.b[i,2])]^0.5)%*%t(phiY.b[,(com.b[i,1])]*lamdaK.b[(com.b[i,1])]^0.5))*Z2.b[b,(com.b[i,1])]
    }

    for(i in 1:kappa.a){
      zeta.a <- zeta.a+(phiY.a[,i]*lamdaK.a[i]^0.5)%*%t(phiY.a[,i]*lamdaK.a[i]^0.5)*((colMeans(epi.a)[i]-1))^0.5*Z1.a[b,i]
    }

    for(i in 1:kappa.b){
      zeta.b <- zeta.b+(phiY.b[,i]*lamdaK.b[i]^0.5)%*%t(phiY.b[,i]*lamdaK.b[i]^0.5)*((colMeans(epi.b)[i]-1))^0.5*Z1.b[b,i]
    }


    Xi[b] <= max(abs(zeta.a-zeta.b)*(V.hat.a - rab*V.hat.b)^(-0.5))
  }
  Q <- quantile(Xi,1-alpha.grid)#quantiles for three levels

  ###   construct the envelope     ###
  ####################################
  sce.band <- array(NA,dim=c(NN, NN,2,nalpha))
  for(j in 1:nalpha){
    sce.band[,,1,j] <- GhatY.a-GhatY.b-Q[j]*na^(-0.5)*(V.hat.a+rab*V.hat.b)^0.5#lower bound
    sce.band[,,2,j] <- GhatY.a-GhatY.b+Q[j]*na^(-0.5)*(V.hat.a+rab*V.hat.b)^0.5#upper bound
  }
  list(sce = sce.band)
}
