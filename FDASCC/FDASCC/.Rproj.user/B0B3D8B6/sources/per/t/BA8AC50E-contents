#' Simultaneous confidence corridors for mean function and its derivatives of one-dimentional functional data
#'
#' The function is used to construct SCC for mean function and its derivatives of one group of images or
#' difference of mean functions between two sets of images.
#'
#' @param Ya a matrix of one-dimensional dense functional data with each row corresponding to one subject.
#' @param Yb an optional matrix containing the second group of functional data. Default is \code{NULL}.
#' When \code{Yb} is \code{NULL}, a one-group SCC is constructed for the mean function or mean function's derivatives of \code{Ya},
#' otherwise, a two-group SCC is constructed for the difference between the mean functions of \code{Yb} and \code{Ya}.
#' @param X an optional vector containing locations of observations.
#' If not specified, the data are treated to be observed on regular grids in \code{(0,1])}.
#' @param X.band an optional vector specifying the locations to construct scc. If not provided, scc will be generated on \code{X}.
#' @param d.est degree of  B-spline for estimating mean function, default is 3 for cubic splines.
#' @param d.band degree of B-spline on each dimension when constructing tensor-product spline smoothing for covariance function, default is 3.
#' @param derivs order of mean function's derivatives to conduct scc. Default is 1.
#' @param nknots.est.a a vector or an integer specifying the potential numbers of knots when estimating mean function with B-spline.
#' Default is \code{seq(2,min(floor(na/4),20),1)}, where \code{na} is the total number of subjects in the first group.
#' If \code{nknots.est.a} has length more than \code{1}, then GCV is used to select the optimal number of knots.
#' @param nknots.cov.a an integer indicating the number of knots in each dimension when using
#' tensor-product splines to estimate the first group's covariance function.
#' @param nknots.est.b,nknots.cov.b optional information about spline estimation for second set of functional data.
#' @param alpha.grid vector of confidence levels. Default is \code{c(0.1,0.05,0.01)}.
#' @param nboot number of bootstrap to simulate quantiles. Default is 1000.
#'
#' @import MASS
#' @importFrom RSpectra eigs
#' @importFrom pracma Real
#' @import splines2
#' @import gtools
#'
#' @export

scc.1D <- function(Ya, Yb=NULL, X=NULL, X.band=NULL, d.est=3, d.cov=3, derivs = 1,
                   nknots.est.a=NULL, nknots.est.b=NULL, nknots.cov.a = 6, nknots.cov.b = NULL,
                   alpha.grid=c(0.1,0.05,0.01), nboot = 1000){
  # Penalty?
  try(if(is.null(Ya)) stop("The response curves are missing."))
  try(if(derivs < 0) stop("The order of derivative should be a positive integer."))
  N <- ncol(Ya)
  if (!is.null(X)){
    try(if(length(X)!=ncol(Ya)) stop("The numbers of points in the response curves and the location matrix are different."))
  }else{
    X <- seq(1/N,1,1/N)
  }
  if (is.null(X.band)){
    X.band <- X
  }

  this.call <- match.call()
  mfit <- list()

  ## One Group SCB
  if(is.null(Yb)){
    out <- scc1g.1D(Y=Ya, X=X, X.band=X.band, d.est=d.est, d.cov=d.cov, derivs=derivs,
                    nknots.est=nknots.est.a, nknots.cov=nknots.cov.a,
                    alpha.grid = alpha.grid, nboot = nboot) #call "scc1g.1D"
    mfit$Yhat <- out$muhat
    mfit$Yhat.pred <- out$muhat.band
    mfit$Yhat.deriv <- out$muhat.deriv
    mfit$Yhat.deriv.pred <- out$muhat.deriv.band
    mfit$scc <- out$scc; mfit$scc.deriv <- out$scc.deriv
    mfit$sce <- out$sce
    mfit$cover.zero <- out$cover.zero;
    mfit$bw <- out$bw; mfit$bw.deriv <- out$bw.deriv
    mfit$knots.est.a <- out$knots.est.a;
    mfit$knots.est.b <- out$knots.est.b;
    mfit$knots.cov.a <- out$knots.cov.a;
    mfit$knots.cov.b <- out$knots.cov.b;
  }

  ## Two Group SCB
  if (!is.null(Yb)){
    if(is.null(nknots.cov.b)){
      nknots.cov.b <- nknots.cov.a
    }
    out <- scc2g.1D(Ya, Yb, X, X.band, d.est, d.cov, derivs, nknots.est.a, nknots.est.b,
                    nknots.cov.a, nknots.cov.b, alpha.grid, nboot)  #call "scc2g.1D"
    mfit$Yhat <- rbind(out$muhat.a, out$muhat.b)
    mfit$Yhat.pred <- out$muhat.band.b - out$muhat.band.a
    mfit$Yhat.deriv <- rbind(out$muhat.deriv.a, out$muhat.deriv.b)
    mfit$Yhat.deriv.pred <- out$muhat.deriv.band.b - out$muhat.deriv.band.a
    mfit$scc <- out$scc; mfit$scc.deriv <- out$scc.deriv
    mfit$sce <- out$sce
    mfit$cover.zero <- out$cover.zero;
    mfit$bw <- out$bw; mfit$bw.deriv <- out$bw.deriv
    mfit$knots.est.a <- out$knots.est.a;
    mfit$knots.est.b <- out$knots.est.b;
    mfit$knots.cov.a <- out$knots.cov.a;
    mfit$knots.cov.b <- out$knots.cov.b;
  }
  mfit$Ya <- Ya; mfit$Yb <- Yb;
  mfit$X <- X; mfit$X.band <- X.band;
  mfit$d.est <- d.est; mfit$d.cov <- d.cov; mfit$derivs <- derivs;
  mfit$alpha <- alpha.grid
  mfit$call <- this.call;
  class(mfit) <- "func"
  return(mfit)
}


