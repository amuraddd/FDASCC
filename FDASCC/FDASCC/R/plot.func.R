#' Plot estimation and sccs of mean function and its derivatives for functional data
#'
#' The plot method for class 'func'. The function plots estimated mean functions (derivatives) and corresponding simultaneous confidence corridors.
#' @param mfit an "func" object returned from either \code{fit.func.1D} or \code{scc.1D}.
#'
#' @export
plot.func <- function(mfit){
  ## Plot 1: Estimated functions/mean function
  X <- mfit$X
  X.band <- mfit$X.band
  for(ii in 1:nrow(mfit$Yhat)){
    Yhat <- mfit$Yhat[ii,]
    plot(X, Yhat, type = 'l', xlab='t', ylab='fitted function')
    readline(prompt="Press [enter] to continue")
  }
  ## Plot 2: SCC
  if(!is.null(mfit$scc)){
    mhat.pred <- mfit$Yhat.pred
    nalpha <- (dim(mfit$scc))[3]
    for(ii in 1:nalpha){
      scc <- mfit$scc[,,ii]
      plot(X.band, mhat.pred, type = 'l',
           ylim = c(min(scc[,1])-1,max(scc[,2])+1))

      title(paste("Lower SCC when alpha =",mfit$alpha[ii]))
      lines(X.band, scc[,1], col='blue')
      lines(X.band, scc[,2], col='blue')
      readline(prompt="Press [enter] to continue")
    }
  }

  ## Plot 3: Estimated functional derivatives
  for(ii in 1:nrow(mfit$Yhat.deriv)){
    Yhat.deriv <- mfit$Yhat.deriv[ii,]
    plot(X, Yhat.deriv, type = 'l', xlab='t', ylab='fitted derivative')
    readline(prompt="Press [enter] to continue")
  }
  ## Plot 4: SCC for functional derivatives
  if(!is.null(mfit$scc.deriv)){
    mhat.deriv.pred <- mfit$Yhat.deriv.pred
    nalpha <- (dim(mfit$scc))[3]
    for(ii in 1:nalpha){
      scc.deriv <- mfit$scc.deriv[,,ii]
      plot(X.band, mhat.deriv.pred, type = 'l',
           ylim = c(min(scc.deriv[,1])-1,max(scc.deriv[,2])+1))

      title(paste("Lower SCC when alpha =",mfit$alpha[ii]))
      lines(X.band, scc.deriv[,1], col='blue')
      lines(X.band, scc.deriv[,2], col='blue')
      readline(prompt="Press [enter] to continue")
    }
  }

}

