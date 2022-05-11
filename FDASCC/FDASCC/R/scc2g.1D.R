#' @import MASS
#' @importFrom RSpectra eigs_sym
#' @importFrom pracma Real
#' @import splines2
#'


scc2g.1D <- function(Ya, Yb, X, X.band, d.est=3, d.cov=3, derivs=1,
                     nknots.est.a, nknots.est.b, nknots.cov.a, nknots.cov.b,
                     alpha.grid=c(0.1,0.05,0.01), nboot = 2000){
  
  #Construct simultaneous confidence bands in two samples case for the difference of true functional curves 
  #and their derivatives(if required) by the estimates of mean functions, covariance function 
  #and quantile. Construct simultaneous confidence envelopes for covariance function.B-spline method is used to estimate the mean function and convariance function, 
  #  Bootstrap is used to generate quantile.
  na <- nrow(Ya)
  nb <- nrow(Yb)
  rab <- na/nb  #r_hat
  N <- ncol(Ya)
  N.band <- length(X.band)
  nalpha <- length(alpha.grid)
  X.sup <- quantile(c(X, X.band), c(0,1))  # support of functions

  if (is.null(nknots.est.a)){
    nknots.est.a <- seq(2,min(floor(na/4),20),1)
    GCV.flag.a <- TRUE
  }else if (length(nknots.est.a) == 1){
    GCV.flag.a <- FALSE
  }else{
    GCV.flag.a <- TRUE
  }
  if (is.null(nknots.est.b)){
    nknots.est.b <- seq(2,min(floor(nb/4),20),1)
    GCV.flag.b <- TRUE
  }else if (length(nknots.est.b) == 1){
    GCV.flag.b <- FALSE
  }else{
    GCV.flag.b <- TRUE
  }

  ## Step 1. Estimate mean function and derivatives
  
  ##Group 1
  Yam <- matrix(apply(Ya,2,mean),nrow=1) #Mean for jth point
  mfit0.a <- fit.mean.1D(Yam, X, X.band, nknots.est.a, d.est, derivs, X.sup, GCV.flag.a)
  #call "fit.mean.1D" estimae mean function
  knots.GCV.a <- unlist(mfit0.a$knots)
  muhat.a <- mfit0.a$Yhat
  muhat.v.a <- mfit0.a$Yhat.deriv
  muhat.band.a <- mfit0.a$Yhat.pred
  muhat.v.band.a <- mfit0.a$Yhat.deriv.pred
  muhat.mtx.a <- matrix(rep(drop(muhat.a),times=na),nrow=na,byrow=TRUE)

  R0.a <- Ya-muhat.mtx.a  #residuals 
  C.a <- crossprod(R0.a)/na
  #C: N*N matrix C[j,k]<-C[k,j]<-mean((Ya[,j]-muhat.a[j])*(Ya[,k]-muhat.a[k]))
  
  
  
  ##Group 2
  Ybm <- matrix(apply(Yb,2,mean),nrow=1)
  mfit0.b <- fit.mean.1D(Ybm, X, X.band, nknots.est.b, d.est, derivs, X.sup, GCV.flag.b)
  #call "fit.mean.1D" estimae mean function
  knots.GCV.b <- unlist(mfit0.b$knots)
  muhat.b <- mfit0.b$Yhat
  muhat.v.b <- mfit0.b$Yhat.deriv
  muhat.band.b <- mfit0.b$Yhat.pred
  muhat.v.band.b <- mfit0.b$Yhat.deriv.pred
  muhat.mtx.b <- matrix(rep(drop(muhat.b),times=nb),nrow=nb,byrow=TRUE)
  
  R0.b <- Yb-muhat.mtx.b  #residuals 
  C.b <- crossprod(R0.b)/nb
  #C: N*N matrix C[j,k]<-C[k,j]<-mean((Yb[,j]-muhat.b[j])*(Yb[,k]-muhat.b[k]))

  
  ## Step 2. Estimate Ghat and Ghat.v
  
  ##Group 1
  cov.fit.a <- fit.cov.1D(C.a, X, X.band, nknots = nknots.cov.a, d.cov, derivs, X.sup)
  #call "fit.cov.1D" estimate covariance function
  GhatY.a <- cov.fit.a$GhatY
  Ghat.a <- cov.fit.a$Ghat
  Ghat.band.a <- cov.fit.a$Ghat.pred
  knots.cov.a <- cov.fit.a$knots

  
  ##Group 2
  cov.fit.b <- fit.cov.1D(C.b, X, X.band, nknots = nknots.cov.b, d.cov, derivs, X.sup)
  #call "fit.cov.1D" estimate covariance function
  GhatY.b <- cov.fit.b$GhatY
  Ghat.b <- cov.fit.b$Ghat
  Ghat.band.b <- cov.fit.b$Ghat.pred
  knots.cov.b <- cov.fit.b$knots

  Uhat <- Ghat.a + rab * Ghat.b
  Uhat.band <- Ghat.band.a + rab * Ghat.band.b  #covariance function for the difference

  ## Step 3. Get eigenvalue and eigenfunctions
  
  ##Group 1
  Geig.a <- eigs_sym(Ghat.a/N,k=10)
  evalues.a <- Real(Geig.a$values)
  Kappa.a <- sum(cumsum(evalues.a)/sum(evalues.a)<0.95)+1 #creterion 0.95
  lambdaK.hat.a <- evalues.a[1:Kappa.a]
  psi.hat.a <- matrix((Geig.a$vectors)[,1:Kappa.a],ncol=Kappa.a) #eigenfunction
  psi.hat.a <- psi.hat.a*sqrt(N)
  phi.hat.a <- psi.hat.a%*%diag(sqrt(lambdaK.hat.a),nrow=Kappa.a) #sqrt{lambdaK}*eigenfunction

  
  ##Group 2
  Geig.b <- eigs_sym(Ghat.b/N,k=10)
  evalues.b <- Real(Geig.b$values)
  Kappa.b <- sum(cumsum(evalues.b)/sum(evalues.b)<0.95)+1 #creterion 0.95
  lambdaK.hat.b <- evalues.b[1:Kappa.b]
  psi.hat.b <- matrix((Geig.b$vectors)[,1:Kappa.b],ncol=Kappa.b) #eigenfunction
  psi.hat.b <- psi.hat.b*sqrt(N)
  phi.hat.b <- psi.hat.b%*%diag(sqrt(lambdaK.hat.b),nrow=Kappa.b) #sqrt{lambdaK}*eigenfunction
  
  
  ##SCE
  cov.sce =  cov2g.sce(na, nb, N, Kappa.a, Kappa.b, Ya, Yb, muhat.a, muhat.b, Ghat.a, Ghat.b,  GhatY.a, GhatY.b, lambdaK.hat.a, lambdaK.hat.b, rab, nboot)
  sce = cov.sce$sce
  
  ## Step 4. Generate Quantile
  Zk.a <- matrix(rnorm(Kappa.a*nboot,0,1),nboot,Kappa.a)
  Zk.b <- matrix(rnorm(Kappa.b*nboot,0,1),nboot,Kappa.b)
  tmp.hat <- diag(1/sqrt(diag(Uhat)))
  zeta.hat <- tmp.hat%*%(phi.hat.a%*%t(Zk.a)-sqrt(rab)*phi.hat.b%*%t(Zk.b))
  Qalpha.hat <- quantile(apply(abs(zeta.hat),2,max),c(1-alpha.grid))

  ## Step 5. Construct SCB
  ME.band <- matrix(rep(Qalpha.hat,each=N.band),nrow=N.band,ncol=nalpha)*
    matrix(rep(sqrt(diag(Uhat.band))/sqrt(na),times=nalpha),nrow=N.band,ncol=nalpha)
  muband.diff <- matrix(rep(muhat.band.b-muhat.band.a,times=nalpha),nrow=N.band,ncol=nalpha) #difference of estimated mean function
  scc.l.band <- muband.diff-ME.band  #lower band
  scc.u.band <- muband.diff+ME.band  #upper band
  bw.band <- 2*apply(ME.band,2,mean)

  #######################################################
  ## SCB for functional derivatives - Two Sample Cases? #
  #######################################################
  if (derivs > 0){
    
    ##Group 1
    Ghat.v.a <- cov.fit.a$Ghat.v
    Ghat.v.band.a <- cov.fit.a$Ghat.v.pred
    psi.v.hat.a <- t(t(psi.hat.a) %*% (Ghat.v.a/N) / lambdaK.hat.a) 
    psi.v.band.a <- t(t(psi.hat.a) %*% (Ghat.v.band.a/N) / lambdaK.hat.a) #eigenfunctions
    phi.v.hat.a <- psi.v.hat.a%*%diag(sqrt(lambdaK.hat.a),nrow=Kappa.a) 
    phi.v.band.a <- psi.v.band.a%*%diag(sqrt(lambdaK.hat.a),nrow=Kappa.a) #eigenfunctions

    
    ##Group 2
    Ghat.v.b <- cov.fit.b$Ghat.v
    Ghat.v.band.b <- cov.fit.b$Ghat.v.pred
    psi.v.hat.b <- t(t(psi.hat.b) %*% (Ghat.v.b/N) / lambdaK.hat.b)  
    psi.v.band.b <- t(t(psi.hat.b) %*% (Ghat.v.band.b/N) / lambdaK.hat.b)  #eigenfunctions
    phi.v.hat.b <- psi.v.hat.b%*%diag(sqrt(lambdaK.hat.b),nrow=Kappa.b)  
    phi.v.band.b <- psi.v.band.b%*%diag(sqrt(lambdaK.hat.b),nrow=Kappa.b)  #eigenfunctions

    Vhat.a <- phi.v.hat.a%*%t(phi.v.hat.a)
    Vhat.band.a <- phi.v.band.a%*%t(phi.v.band.a) #sigma function 

    Vhat.b <- phi.v.hat.b%*%t(phi.v.hat.b)
    Vhat.band.b <- phi.v.band.b%*%t(phi.v.band.b) #sigma function 

    What <- Vhat.a + rab * Vhat.b
    What.band <- Vhat.band.a + rab * Vhat.band.b #sigma function for the difference

    tmp.v.hat <- diag(1/sqrt(diag(What)))
    zeta.v.hat <- phi.v.hat.a%*%t(Zk.a)-sqrt(rab)*phi.v.hat.b%*%t(Zk.b)
    zeta.v.hat <- tmp.v.hat%*%zeta.v.hat
    Qalpha.v.hat <- quantile(apply(abs(zeta.v.hat),2,max),c(1-alpha.grid)) #generated quantile

    ME.v.band <- matrix(rep(Qalpha.v.hat,each=N.band),nrow=N.band,ncol=nalpha)*
      matrix(rep(sqrt(diag(What.band))/sqrt(na),times=nalpha),nrow=N.band,ncol=nalpha)
    muband.v.diff <- matrix(rep(muhat.v.band.b-muhat.v.band.a,times=nalpha),nrow=N.band,ncol=nalpha)
    scc.l.v.band <- muband.v.diff-ME.v.band  #lower band
    scc.u.v.band <- muband.v.diff+ME.v.band  #upper band
    bw.v.band <- 2*apply(ME.v.band,2,mean)
  }else{
    scc.l.v.band <- scc.l.band
    scc.u.v.band <- scc.u.band
    bw.v.band <- bw.band
  }
  scc.band <- array(NA,dim=c(N.band,2,nalpha))
  scc.v.band <- array(NA,dim=c(N.band,2,nalpha))
  cover.zero <- matrix(NA,nrow=N.band,ncol=nalpha)
  cover.zero.v <- matrix(NA,nrow=N.band,ncol=nalpha)
  for(j in 1:nalpha){
    scc.band[,,j] <- cbind(scc.l.band[,j],scc.u.band[,j])
    scc.v.band[,,j] <- cbind(scc.l.v.band[,j],scc.u.v.band[,j])
    cover.zero[,j] <- matrix(as.numeric((0>scc.u.band[,j]))-as.numeric((0<scc.l.band[,j])))
    cover.zero.v[,j] <- matrix(as.numeric((0>scc.u.v.band[,j]))-as.numeric((0<scc.l.v.band[,j])))
  }
  list(nknots.est.a = mfit0.a$Nsc, nknots.est.b = mfit0.b$Nsc, knots.est.a = knots.GCV.a, knots.est.b = knots.GCV.b,
       knots.cov.a = knots.cov.a, knots.cov.b = knots.cov.b,
       d.est = d.est, d.cov = d.cov,
       muhat.a = muhat.a, muhat.b = muhat.b, muhat.band.a = muhat.band.a, muhat.band.b = muhat.band.b,
       muhat.deriv.a = muhat.v.a, muhat.deriv.b = muhat.v.b,
       muhat.deriv.band.a = muhat.v.band.a, muhat.deriv.band.b = muhat.v.band.b,
       scc = scc.band, scc.deriv = scc.v.band, sce = sce,
       cover.zero = cover.zero, cover.zero.deriv = cover.zero.v,
       Ghat.a = Ghat.a, Ghat.b = Ghat.b, Vhat.a = Vhat.a, Vhat.b = Vhat.b,
       bw = bw.band, bw.deriv = bw.v.band)
}


