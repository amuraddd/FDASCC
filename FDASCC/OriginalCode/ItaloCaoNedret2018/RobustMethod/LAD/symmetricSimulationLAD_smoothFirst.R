########
#Confidence Band Test - Spline band----
########


source("auxCode/gervini_functions.R")
source("auxCode/auxFunctions.R")
#Libraries ----
library(splines)
library(MASS)
library(quantreg)
library(reshape2)
#Sample Generation ----

symmetricSimLADSmoothFirst <- function(contamination, n, p = 4, Nm = 0, c = 0.5, sigma = 0.5, alpha = 0.05, strength = 10, correction = 0, debug = FALSE) { 
  #Constants ----
  
  # contamination = Proportion of dataset of outliers
  # n             = Number of samples
  # p - 1         = Degree of spline polynomials
  # Nm            = Number of knots (if 0, use formula from Cao)
  # c             = Constant used to find Nm
  # sigma         = Variance of the measrument error
  # alpha         = Confidence band
  # strength      = strengh of the outlier - dependent on type of outlier
  
  N <- floor(n/2) #Number of observations
  
  if(Nm == 0) {
    Nm <- floor(c*(n^(1/(2*p+1))*log(n))) #Number of  knots
  }
  
  #Definition of Correction Factor
  if(correction == 0) {
    # correction = 1.3*1.1^(log(n/30, base = 2))
    correction = 3/(1+n^0.05)
  }
  #####
  #####
  
  N.test = 100
  
  #Spline pre-calculations ----
  x <- (1:N)/N
  knots<- generateKnots(p, Nm, x)
  
  #Pre-calculating the B-Splines 
  
  x.melt <- rep(x,n)
  
  X.long <- bs(x = x.melt , knots =  knots[2:(Nm+1)], degree = (p-1), intercept=TRUE, Boundary.knots=knots[c(1,Nm+2)])
  X <- bs(x = x , knots =  knots[2:(Nm+1)], degree = (p-1), intercept=TRUE, Boundary.knots=knots[c(1,Nm+2)])
  
  X.pre.calc <-  solve(t(X) %*% X) %*% t(X)
  t.coverage.test <- (1:N.test)/N.test
  
  X.coverage.test <- bs(x = t.coverage.test , knots =  knots[2:(Nm+1)], degree = (p-1), intercept=TRUE, Boundary.knots=knots[c(1,Nm+2)])
  
  real.mean <- generate.mean(N)
  mean.coverage <- generate.mean(N.test)
  
  #Beggining of simulation ----
  
  miss <- 0
  
  n.out = floor(n*contamination)
  n.real = n - n.out
  maxSim <- 100
  
  if (debug == TRUE) {
    l2mean <- 0
    bandSize <- 0 
  }
  pb <- txtProgressBar(min = 0, max = maxSim)
  
  for (a in 1 : maxSim){
    
    sample <- rbind(generateSample1(n.real, N, sigma), generateOutlierSymmetric(numberOfSamples =  n.out, N, strength = strength, sigma))
    sample.melt <- melt(t(sample))$value
    
    #Mean function Estimation ----  
    coef.mean <- as.numeric(rq(sample.melt ~ X.long - 1)$coef)
    mean.hat<- X %*% coef.mean
    
    #Covariance Estimation ----  
    sp.princomp <- f.Sp.princomp.2(t = x, x = sample, mu = as.vector(mean.hat))
    
    #Eigenvalue calculation in R is resturning negative values, instead of ~ 0, in the middle of algorithm- This amounts to returning NA for certaing scores
    # sp.princomp$scores[is.na(sp.princomp$scores)] = 0
    # eigenvalues <- apply(sp.princomp$scores, 2, function(data) hubers(data, mu = 0)$s)^2
    
    eigenvectors <- sp.princomp$pc
    eigenvectors[is.na(eigenvectors)] = 0
    eigenvectors.coef <- X.pre.calc %*% t(eigenvectors)
    eigenvectors <- t(X %*% eigenvectors.coef)
    sample.centered <- t(t(sample) - as.vector(mean.hat))
    eigenvalues <- apply(sample.centered %*% t(eigenvectors)/N, 2, function(data) hubers(data , mu=0)$s)^2
    
    eigenvalues.sum <- sum(eigenvalues, na.rm = T)
    
    for (i in 1:length(eigenvalues)){
      #Find the smallest number of eigenvalues that will make up for 95% of the total
      if(sum(eigenvalues[1:i])/eigenvalues.sum > 0.95)
      {
        k <- i
        break
      }
    }
    

    eigenvalues <- eigenvalues[1:k]
    if(k == 1){
      eigenvectors <- t(as.matrix(eigenvectors[1, ]))
    } else {
      eigenvectors <-eigenvectors[1:k , ]
    }    
    eigenvectors.coef <- X.pre.calc %*% t(eigenvectors)
    eigenvectors <- t(X %*% eigenvectors.coef)
    G.hat <- t(eigenvectors) %*%(eigenvalues*eigenvectors)
    
    
    phi.hat = matrix()
    length(phi.hat) <- k*N
    dim(phi.hat) <- c(k,N)
    for (i in 1:k){
      phi.hat[i,] <- sqrt(eigenvalues[i])*eigenvectors[i,]
    }  
    
    
    
    max.sim <- array(data=0.0, dim=1000)
    for(i in 1: 1000){
      z = rnorm(dim(phi.hat)[1])
      max.sim[i] = max(colSums(z*phi.hat)/sqrt(diag(abs(G.hat))))
    }
    
    Q.hat = quantile(max.sim, probs = c(1-alpha))
    
    mean.hat.coverage.test <- X.coverage.test%*%coef.mean 
    
    eigenvectors.test <- t(X.coverage.test%*%eigenvectors.coef)
    
    G.hat.test <- t(eigenvectors.test) %*%(eigenvalues*eigenvectors.test)
    dim(G.hat.test) <- c(N.test, N.test)
    
    band <- correction*(sqrt(diag(abs(G.hat.test)))*Q.hat)/sqrt(n)
    
    upper.limit = mean.hat.coverage.test + band
    lower.limit = mean.hat.coverage.test - band
    #Coverage calculation ----
    
    if(all(mean.coverage <= upper.limit) & all(mean.coverage >= lower.limit))
      miss <- miss + 1
    
    if(debug == TRUE) {
      l2mean <- l2mean + abs(mean.hat.coverage.test - mean.coverage)
      bandSize <-  bandSize + abs(band)
    }
    setTxtProgressBar(pb, a)
    
  }
  close(pb)
  #   cat( proc.time() - timeBefore )
  if(debug == TRUE){
    return(list(l2mean = l2mean/maxSim, bandSize = bandSize/maxSim, coverage = miss/maxSim))
  } else {
    return(miss/maxSim)  
  }
  
  
  
}

