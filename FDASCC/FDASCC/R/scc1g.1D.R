#' @import MASS
#' @importFrom RSpectra eigs_sym
#' @importFrom pracma Real
#' @import splines2
#'


scc1g.1D <- function(Y, X, X.band, d.est, d.cov, derivs=1, nknots.est, nknots.cov=6,
                   alpha.grid=c(0.1,0.05,0.01), nboot = 2000){
  #Construct simultaneous confidence bands for the both true functional curve and 
  #its derivatives(if required) by the estimates of mean function, covariance function 
  #and quantile. Construct simultaneous confidence envelopes for covariance function. B-spline is used to estimate the mean function and convariance function, 
  #  Bootstrap is used to generate quantile.
  n <- nrow(Y)
  N <- ncol(Y)
  N.band <- length(X.band)
  nalpha <- length(alpha.grid)
  X.sup <- quantile(c(X, X.band), c(0,1))  # support of functions

  if (is.null(nknots.est)){
    nknots.est <- seq(2,min(floor(n/4),20),1)
    GCV.flag <- TRUE
  }else if (length(nknots.est) == 1){
    # Ns1 <- nknots.est
    GCV.flag <- FALSE
  }else{
    GCV.flag <- TRUE
  }

  ## Step 1. Estimate mean function and derivatives
  Ym <- matrix(apply(Y,2,mean),nrow=1)
  mfit0 <- fit.mean.1D(Ym, X, X.band, nknots=nknots.est, d=d.est, derivs, X.sup, GCV.flag)
  #call "fit.mean.1D" estimae mean function
  knots.GCV <- unlist(mfit0$knots)
  muhat <- mfit0$Yhat
  muhat.v <- mfit0$Yhat.deriv
  muhat.band <- mfit0$Yhat.pred
  muhat.v.band <- mfit0$Yhat.deriv.pred
  muhat.mtx <- matrix(rep(drop(muhat),times=n),nrow=n,byrow=TRUE)

  R0 <- Y-muhat.mtx #residuals
  C <- crossprod(R0)/n
  #C: N*N matrix C[j,k]<-C[k,j]<-mean((Y[,j]-muhat[j])*(Y[,k]-muhat[k]))

  ## Step 2. Estimate Ghat and Ghat.v
  cov.fit <- fit.cov.1D(C, X, X.band, nknots = nknots.cov, d.cov, derivs, X.sup)
  GhatY <- cov.fit$GhatY
  Ghat <- cov.fit$Ghat
  Ghat.band <- cov.fit$Ghat.pred
  knots.cov <- cov.fit$knots

  ## Step 3. Get eigenvalue and eigenfunctions
  Geig <- eigs_sym(Ghat/N,k=10)
  evalues <- Real(Geig$values)
  Kappa <- sum(cumsum(evalues)/sum(evalues)<0.95)+1  #creterion 0.95
  lambdaK.hat <- evalues[1:Kappa]
  psi.hat <- matrix((Geig$vectors)[,1:Kappa],ncol=Kappa) #eigenfunction
  psi.hat <- psi.hat*sqrt(N)
  phi.hat <- psi.hat%*%diag(sqrt(lambdaK.hat),nrow=Kappa) #sqrt{lambdaK}*eigenfunction

  ## SCE
 
  cov.sce = cov1g.sce(n, N, Kappa, Y, muhat, Ghat, GhatY, lambdaK.hat, nboot)
  sce = cov.sce$sce
  
  ## Step 4. Generate Quantile for SCB
  Zb <- matrix(rnorm(Kappa*nboot,0,1),nboot,Kappa)
  tmp.hat <- diag(1/sqrt(diag(Ghat)))
  zeta.hat <- phi.hat%*%t(Zb)
  zeta.hat <- tmp.hat%*%zeta.hat
  Qalpha.hat <- quantile(apply(abs(zeta.hat),2,max),c(1-alpha.grid))

  ## Step 5. Construct SCB
  ME.band <- matrix(rep(Qalpha.hat,each=N.band),nrow=N.band,ncol=nalpha)*
    matrix(rep(sqrt(diag(Ghat.band))/sqrt(n),times=nalpha),nrow=N.band,ncol=nalpha)
  muband.mtx <- matrix(rep(muhat.band,times=nalpha),nrow=N.band,ncol=nalpha)
  scc.l.band <- muband.mtx-ME.band #lower band
  scc.u.band <- muband.mtx+ME.band #upper band
  bw.band <- 2*apply(ME.band,2,mean)

  ####################################
  ## SCB for functional derivatives
  ####################################
  if (derivs > 0){
    Ghat.v <- cov.fit$Ghat.v
    Ghat.v.band <- cov.fit$Ghat.v.pred
    psi.v.hat <- t(t(psi.hat) %*% (Ghat.v/N) / lambdaK.hat)
    psi.v.band <- t(t(psi.hat) %*% (Ghat.v.band/N) / lambdaK.hat)
    phi.v.hat <- psi.v.hat%*%diag(sqrt(lambdaK.hat),nrow=Kappa)
    phi.v.band <- psi.v.band%*%diag(sqrt(lambdaK.hat),nrow=Kappa)
    
    Vhat <- phi.v.hat%*%t(phi.v.hat)
    Vhat.band <- phi.v.band%*%t(phi.v.band) #sigma function

    tmp.v.hat <- diag(1/sqrt(diag(Vhat)))
    zeta.v.hat <- phi.v.hat%*%t(Zb)
    zeta.v.hat <- tmp.v.hat%*%zeta.v.hat
    Qalpha.v.hat <- quantile(apply(abs(zeta.v.hat),2,max),c(1-alpha.grid)) #generated quantile

    ME.v.band <- matrix(rep(Qalpha.v.hat,each=N.band),nrow=N.band,ncol=nalpha)*
      matrix(rep(sqrt(diag(Vhat.band))/sqrt(n),times=nalpha),nrow=N.band,ncol=nalpha)
    muband.v.mtx <- matrix(rep(muhat.v.band,times=nalpha),nrow=N.band,ncol=nalpha)
    scc.l.v.band <- muband.v.mtx-ME.v.band #lower band
    scc.u.v.band <- muband.v.mtx+ME.v.band #upper band
    bw.v.band <- 2*apply(ME.v.band,2,mean)
  }else{
    scc.l.v.band <- scc.l.band
    scc.u.v.band <- scc.u.band
    bw.v.band <- bw.band
  }
  scc.band <- array(NA,dim=c(N.band,2,nalpha))
  scc.v.band <- array(NA,dim=c(N.band,2,nalpha))
  for(j in 1:nalpha){
    scc.band[,,j] <- cbind(scc.l.band[,j],scc.u.band[,j])
    scc.v.band[,,j] <- cbind(scc.l.v.band[,j],scc.u.v.band[,j])
  }
  list(nknots.est = mfit0$Nsc, knots.est = knots.GCV, knots.cov = knots.cov,
       d.est = d.est, d.cov = d.cov,
       muhat = muhat, muhat.band = muhat.band,
       muhat.deriv = muhat.v, muhat.deriv.band = muhat.v.band,
       scc = scc.band, scc.deriv = scc.v.band,
       Ghat = Ghat, Vhat = Vhat, bw = bw.band, bw.deriv = bw.v.band, sce = sce)
}




