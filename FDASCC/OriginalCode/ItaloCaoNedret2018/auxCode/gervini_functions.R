

#Functional median
f.Sp.median <- function(t, x, bias.corr = FALSE) {
  #x - Sample in each line - n,m
  #t - time points - m
  
  n <- nrow(x)
  m <- ncol(x)
  
  A  = (x[ , 1:(m-1)] %*% diag(t[2:m] - t[1:(m-1)]) %*% t(x[ , 1:(m-1)]) + 
          x[ , 2:m] %*% diag(t[2:m] - t[1:(m-1)]) %*% t(x[ , 2:m]))*0.5
  
  if(bias.corr) {
    A = A - diag(gasser(t, x), nrow = n)*(t[m] - t[1])
  }
  
  w <- rep(1/n, n)
  norms <- sqrt( diag(A) + t(w) %*% A %*% w - 2* A %*% w)
  f <- sum(norms)
  err <- 1
  iter <- 0
  
  while( err > 1e-5 && iter < 50) {
    iter <- iter + 1
    f0 <- f
    if (any(norms < .Machine$double.eps)) {
      i0 <- norms < .Machine$double.eps
      w <- rep(0, n)
      w[i0] <- 1/length(i0)
    } else {
      w <- 1/norms
      w <- w/sum(w)
    }
    
    norms <- sqrt( diag(A) + t(w) %*% A %*% w - 2* A %*% w)
    f <- sum(norms)
    err <- abs(f/f0 - 1)
  }
  
  list(median = as.vector(t(w) %*% x), weights =  w, norms =  norms);
}

# Uses results from functional median
f.Sp.princomp <- function(t, x, wmu, bias.corr = FALSE) {
#   require(rARPACK)
  n <- nrow(x)
  m <- ncol(x)
  
  A  = (x[ , 1:(m-1)] %*% diag(t[2:m] - t[1:(m-1)]) %*% t(x[ , 1:(m-1)]) + 
          x[ , 2:m] %*% diag(t[2:m] - t[1:(m-1)]) %*% t(x[ , 2:m]))*0.5
  
  if(bias.corr) {
    A = A - diag(gasser(t, x), nrow = n)*(t[m] - t[1])
  }
  
  B <- ( diag(1, n, n) - matrix(1, n , 1) %*% t(wmu) ) %*% A %*% (  diag(1, n, n) - wmu %*% matrix(1, 1, n) )
  
  norms <- sqrt(diag(B))
  
  DN <- rep(0, n)
  
  I <- which(norms > .Machine$double.eps)
  DN[I] <- 1 / norms[I]
  B <- diag(DN) %*% B %*% diag(DN)
  B <- (B + t(B))/2
  
  
  eig.decomp <- eigen(B, symmetric = TRUE)
  D <- abs(eig.decomp$values)
  V <- eig.decomp$vectors

#   eig.decomp = eigs_sym(A = B, k = q,which =  "LA")
#   D = eig.decomp$values
#   V = eig.decomp$vectors
#   
  #   sorting.res <- sort(D, decreasing = TRUE, index.return = TRUE)
  #   D <- sorting.res$x
  #   indices <- sorting.res$ix
  
  wpc <- V %*% diag(1/sqrt(D))
  lambda <- D/n
  
  pc <- t(wpc) %*% diag(DN) %*% (diag(1, n, n) - rep(1,n) %*% t(wmu)) %*% x
  
  scores <- 0.5 * (pc[ , 1:(m-1)] %*% diag(t[2:m] - t[1:(m-1)]) %*% t((diag(1, n, n) - rep(1,n) %*% t(wmu)) %*% x[ , 1:(m-1)] ) + pc[ , 2:m] %*% diag(t[2:m] - t[1:(m-1)]) %*%  t((diag(1, n, n) - rep(1,n) %*% t(wmu)) %*% x[ , 2:m]) )
  
  scores <- t(scores)
  
  return.result <- list()
  return.result$pc <- pc
  return.result$wpc <- wpc
  return.result$lambda <- lambda
  return.result$scores <- scores
  
  return(return.result)
}

#Uses general mean function - input: mu
f.Sp.princomp.2 <- function(t, x, mu, bias.corr = FALSE) {
  
  n <- nrow(x)
  m <- ncol(x)
  
  y = t(t(x) - mu)
  G  = (y[ , 1:(m-1)] %*% diag(t[2:m] - t[1:(m-1)]) %*% t(y[ , 1:(m-1)]) + 
          y[ , 2:m] %*% diag(t[2:m] - t[1:(m-1)]) %*% t(y[ , 2:m]))*0.5
  

#   
#   if(bias.corr) {
#     A = A - diag(gasser(t, x), nrow = n)*(t[m] - t[1])
#   }
  
  norms <- sqrt(diag(G))
  
  DN <- rep(0, n)
  
  I <- which(norms > .Machine$double.eps)
  DN[I] <- 1 / norms[I]
  G <- diag(DN) %*% G %*% diag(DN)
  G <- (G + t(G))/2
  
  
  eig.decomp <- eigen(G, symmetric = TRUE)
  D <- abs(eig.decomp$values)
  V <- eig.decomp$vectors
  
  #   eig.decomp = eigs_sym(A = B, k = q,which =  "LA")
  #   D = eig.decomp$values
  #   V = eig.decomp$vectors
  #   
  #   sorting.res <- sort(D, decreasing = TRUE, index.return = TRUE)
  #   D <- sorting.res$x
  #   indices <- sorting.res$ix
  
  wpc <- V %*% diag(1/sqrt(D))
  lambda <- D/n
  
  pc <- t(wpc) %*% diag(DN) %*% y
  
  scores <- 0.5 * (pc[ , 1:(m-1)] %*% diag(t[2:m] - t[1:(m-1)]) %*% t(y[ , 1:(m-1)] ) + pc[ , 2:m] %*% diag(t[2:m] - t[1:(m-1)]) %*%  t(y[ , 2:m]) )
  
  scores <- t(scores)
  
  return.result <- list()
  return.result$pc <- pc
  return.result$wpc <- wpc
  return.result$lambda <- lambda
  return.result$scores <- scores
  
  return(return.result)
}

# Not currently used ----
####
gasser = function(t, x) {
  N = length(t)
  a = (t[3:N] - t[2:(N-1)])/(t[3:N] - t[1:(N-2)])
  b = (t[2:(N-1)] - t[1:(N-2)])/(t[3:N] - t[1:(N-2)])
  c = (a^2 + b^2 + 1)
  
  #return
  median(rowMeans(t(a*t(x[, 1:(N-2)]) + b*t(x[ ,3:N]) - t(x[ ,2:(N-1)]))^2/c))
}
