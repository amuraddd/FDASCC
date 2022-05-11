#' Estimate the functions and functional derivatives via B-spline
#'
#' The function is used to estimate mean function of one-dimensional functional data and its derivatives.
#' The number of knots is selected using generalized cross validation.
#'
#' @param Y a matrix of one-dimensional dense functional data with each row corresponding to one subject.
#' @param X an optional vector containing locations of observations.
#' If not specified, the data are treated to be observed on regular grids in \code{(0,1])}.
#' @param X.pred an optional vector specifying the locations on which the estimated mean or derivative
#' functions are returned in the output. If not provided, estimations will returned on the same locations as in \code{X}.
#' @param nknots a vector or an integer specifying the potential numbers of knots used in B-spline estimation.
#' Default is \code{seq(2,min(floor(n/4),20),1)}, where \code{n} is the total number of subjects.
#' @param d degree of  B-spline for estimating mean function, default is 3 for cubic splines.
#' @param derivs order of functional derivatives to estimate. Default is 1.
#' @param X.sup a length 2 vector of the boundary knots. Default is minimum and maximum in \code{c(X, X.pred)}.
#' @param GCV.flag a logical value indicating whether generalized cross validation should be used to select \code{nknots}.
#'
#' @import splines2
#'
#' @export

fit.func.1D <- function(Y, X, X.pred=NULL, nknots=NULL, d=3, derivs=1, X.sup=NULL, GCV.flag=TRUE){
  # Y is a n*N matrix
  if(!is.matrix(Y)){
    warning("The response variable, Y, should be a matrix with each row represents an image.")
    Y <- as.matrix(Y)
  }
  n <- nrow(Y)
  N <- ncol(Y)
  if (is.null(nknots)){
    nknots <- seq(2,min(floor(n/4),20),1)
  }
  if (is.null(X.pred)){
    X.pred <- X
  }
  if (is.null(X.sup)){
    X.sup <- quantile(c(X, X.pred), c(0,1))
  }
  if (GCV.flag){
    GCV.all <- sapply(nknots, FUN = function(Ns){
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
    Ns.all <- nknots[apply(GCV.all, 1, which.min)]
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
  mfit <- list()
  mfit$theta <- theta;
  mfit$Yhat <- Yhat
  mfit$Yhat.pred <- Yhat.pred
  mfit$Yhat.deriv <- Yhat.v
  mfit$Yhat.deriv.pred <- Yhat.v.pred
  mfit$scc <- NULL;
  mfit$scc.deriv <- NULL;
  mfit$X <- X; mfit$X.band <- X.pred;
  mfit$d.est <- d.est; mfit$d.cov <- d.cov; mfit$derivs <- derivs;
  mfit$alpha <- alpha.grid
  mfit$call <- this.call;
  class(mfit) <- "func"
  return(mfit)
}
