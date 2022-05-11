#' @import splines2
#'
#' @export

fit.cov.1D <- function(C, X, X.pred=NULL, nknots, d, derivs, X.sup){
  #Estimate covariance function for both functional and its derivatives
  #Specific points can be chosen for prediction.
  
  # C: N*N matrix C[j,k]<-C[k,j]<-mean((YY[,j]-mhat[j])*(YY[,k]-mhat[k]))
  
  
  NN=200  #test points for SCE
  ##Part I: Estimation
  
  N <- length(X)
  remov <- seq(1,(N*(N-1)+1), by=N)
  #index for points to be removed
  Cij <- matrix(C,N^2,1)
  Cij <- Cij[-remov,]
  #remove diagonal elements
  Ns <- nknots
  knots <- seq(X.sup[1],X.sup[2],length.out=Ns+2)
  XB <- bSpline(X, knots = knots[-c(1,(Ns+2))], degree = d, intercept = TRUE,
                Boundary.knots = knots[c(1,(Ns+2))])
  TP.XB <- kronecker(XB,XB)
  #B-splines tensor product 
  TP.XB.r <- TP.XB[-remov,]
  #revmoe diagonal basis
  Beta <- solve(t(TP.XB.r)%*%(TP.XB.r))%*%t(TP.XB.r)
  theta <- Beta %*% Cij
  Ghatij <- TP.XB %*% theta
  Ghat <- matrix(Ghatij, N, N, byrow = TRUE)
  #Test points for SCE
  test <- seq(X.sup[1],X.sup[2],length.out = NN)
  XB.test = bSpline(test, knots = knots[-c(1,(Ns+2))], degree = d, intercept = TRUE,
                    Boundary.knots = knots[c(1,(Ns+2))])
  TP.XB.test <- kronecker(XB.test,XB.test)
  Ghatij.test <- TP.XB.test %*% theta
  GhatY <- matrix(Ghatij.test, NN, NN, byrow = TRUE)
  #estimate of covariance function

  if (derivs > 0){
    XB.v <- dbs(X, derivs = derivs, knots = knots[-c(1,(Ns+2))], degree = d,
                    intercept = TRUE, Boundary.knots = knots[c(1,(Ns+2))])
    TP.XB.v <- kronecker(XB,XB.v)
    #B-splines tensor product for derivatives
    Ghatij.v <- TP.XB.v %*% theta
    Ghat.v <- matrix(Ghatij.v, N, N, byrow = TRUE)
    #estimate of covariance function for derivatives
  }else if (derivs == 0){
    Ghat.v <- Ghat
  }

  
  
  
  ##Part II: Prediction
  if (!is.null(X.pred)){
    N.pred <- length(X.pred)
    XB.pred <- bSpline(X.pred, knots = knots[-c(1,(Ns+2))], degree = d, intercept = TRUE,
                  Boundary.knots = knots[c(1,(Ns+2))])
    TP.XB.pred <- kronecker(XB.pred,XB.pred)
    Ghatij.pred <- TP.XB.pred %*% theta
    Ghat.pred <- matrix(Ghatij.pred, N.pred, N.pred, byrow = TRUE)
    #Predicted convariance at X.pred

    if (derivs > 0){
      XB.v.pred <- dbs(X.pred, derivs = derivs, knots = knots[-c(1,(Ns+2))], degree = d,
                  intercept = TRUE, Boundary.knots = knots[c(1,(Ns+2))])
      TP.XB.v.pred <- kronecker(XB,XB.v.pred)
      Ghatij.v.pred <- TP.XB.v.pred %*% theta
      Ghat.v.pred <- matrix(Ghatij.v.pred, N, N.pred, byrow = TRUE)
      #Predicted convariance for derivatives at X.pred 
    }else if (derivs == 0){
      Ghat.v.pred <- Ghat.pred
    }
  }else{
    X.pred <- X
    Ghat.pred <- Ghat
    Ghat.v.pred <- Ghat.v
  }
  list(theta = theta, knots = knots, Ghat = Ghat, GhatY = GhatY, Ghat.v = Ghat.v,
       Ghat.pred = Ghat.pred, Ghat.v.pred = Ghat.v.pred)
}




