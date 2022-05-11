#' @import splines2
#'
#' @export

fit.mean.1D <- function(Y, X, X.pred=NULL, nknots, d, derivs, X.sup, GCV.flag){
  #Estimate the mean function of one-dimensional and its derivatives, B-splines are used as basis functions
  #GCV is used to select knots.
  
  
  # Y is a n*N matrix
  if(!is.matrix(Y)){
    warning("The response variable, Y, should be a matrix with each row represents an image.")
    Y <- as.matrix(Y)
  }
  n <- nrow(Y)
  N <- ncol(Y)
  if (is.null(X.pred)){
    #if there are some specific points to predict
    X.pred <- X
  }
  if (GCV.flag){
    GCV.all <- sapply(nknots, FUN = function(Ns){
      #GCV selection for optimal knots
      knots <- seq(X.sup[1],X.sup[2],length.out=Ns+2)
      XB <- bSpline(X, knots = knots[-c(1,(Ns+2))], degree = d, intercept = TRUE, Boundary.knots = knots[c(1,(Ns+2))])
      beta <- solve(t(XB)%*%(XB))%*%t(XB)
      Hlamda <- XB %*% beta
      df <- sum(diag(Hlamda))
      Yhat <- Hlamda %*% t(Y)
      res <- t(Y)-Yhat
      sse <- apply(res^2,2,sum)   # sse w.r.t. n curves
      gcv <- N*sse/(N-df)^2   
      return(gcv)
    })   # A n * length(nknots) matrix
    GCV.all <- matrix(GCV.all, nrow = n)
    Ns.all <- nknots[apply(GCV.all, 1, which.min)] # Optimal knots number w.r.t n curves
  }else{
    Ns.all <- nknots[1]
    if (length(nknots)>1) {
      warning('The first element of nknots is directly used when GCV.flag==FALSE')
    }
  }

  # Use selected knots
  Yhat <- c()
  Yhat.pred <- c()
  Yhat.v <- c()
  Yhat.v.pred <- c()
  df <- c()
  knots <- vector('list', length = n)
  theta <- vector('list', length = n)
  for (i in 1:n){
    #Prediction
    Ns.tmp <- Ns.all[i]
    knots.tmp <- seq(X.sup[1],X.sup[2],length=Ns.tmp+2)
    knots[[i]] <- knots.tmp
    XB.tmp <- bSpline(X, knots = knots.tmp[-c(1,(Ns.tmp+2))], degree = d, intercept = TRUE,
                      Boundary.knots = knots.tmp[c(1,(Ns.tmp+2))])
    beta.tmp <- solve(t(XB.tmp)%*%(XB.tmp))%*%t(XB.tmp)
    Hlamda.tmp <- XB.tmp %*% beta.tmp
    df.tmp <- sum(diag(Hlamda.tmp))
    df <- c(df, df.tmp)
    Yi <- as.matrix(Y[i,], ncol = 1)
    theta.tmp <- beta.tmp %*% Yi
    theta[[i]] <- theta.tmp
    Yihat <- XB.tmp %*% theta.tmp
    Yhat <- cbind(Yhat, Yihat)
    XB.tmp.pred <- bSpline(X.pred, knots = knots.tmp[-c(1,(Ns.tmp+2))], degree = d, intercept = TRUE,
                      Boundary.knots = knots.tmp[c(1,(Ns.tmp+2))])
    Yihat.pred <- XB.tmp.pred %*% theta.tmp
    Yhat.pred <- cbind(Yhat.pred, Yihat.pred)
    if (derivs > 0){
      #Prediction for derivatives
      XB.v.tmp <- dbs(X, derivs = derivs, knots = knots.tmp[-c(1,(Ns.tmp+2))], degree = d,
                      intercept = TRUE, Boundary.knots = knots.tmp[c(1,(Ns.tmp+2))])
      Yihat.v <- XB.v.tmp %*% theta.tmp
      Yhat.v <- cbind(Yhat.v, Yihat.v)
      XB.v.tmp.pred <- dbs(X.pred, derivs = derivs, knots = knots.tmp[-c(1,(Ns.tmp+2))], degree = d,
                      intercept = TRUE, Boundary.knots = knots.tmp[c(1,(Ns.tmp+2))])
      Yihat.v.pred <- XB.v.tmp.pred %*% theta.tmp
      Yhat.v.pred <- cbind(Yhat.v.pred, Yihat.v.pred)
    }
  }
  Yhat <- t(Yhat)
  Yhat.pred <- t(Yhat.pred)
  Yhat.v <- t(Yhat.v)
  Yhat.v.pred <- t(Yhat.v.pred)
  res <- Y-Yhat
  sse <- apply(res^2,1,sum)
  GCV <- N*sse/(N-df)^2

  list(theta = theta, Yhat = Yhat, Yhat.deriv = Yhat.v,
       Yhat.pred = Yhat.pred, Yhat.deriv.pred = Yhat.v.pred,
       sse = sse, GCV = GCV, Nsc = Ns.all, knots = knots)
}
